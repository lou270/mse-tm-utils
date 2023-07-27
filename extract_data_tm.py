import csv
import matplotlib.pyplot as plt
import numpy as np
import simplekml
import math
from time import time
from datetime import datetime
from scipy.interpolate import CubicSpline

### Extract CSV data
raw_data = []
with open("data/data_tm_simplified.csv", 'r') as fp:
    csvreader = csv.reader(fp, delimiter=';')

    firstLineFound = False
    for row in csvreader:
        if (firstLineFound == False):
            field_name = row
            firstLineFound = True
        else:
            raw_data.append(row)

### Usefuls variables
const_rocket_states = ["PRE_FLIGHT", "PYRO_RDY", "ASCEND", "DEPLOY_ALGO", "DEPLOY_TIMER", "DESCEND", "TOUCHDOWN"]
const_gnss_valid = ["INVALID", "VALID"]
acc_factor = 16.0/32768.0
pitot_factor = 1

# Pitot calibration
calib_pression_pitot = np.array([1.00, 1.04, 1.07, 1.11, 1.16, 1.20, 1.25, 1.31, 1.37, 1.43, 1.50, 1.58, 1.67, 1.77, 1.88, 2.00, 2.15, 2.31, 2.50, 2.73, 3.00, 3.34, 3.76, 4.29, 5.01, 6.01, 7.51]) * 1000 - (1146-968)
calib_adc_pitot = np.array([-24, -30, -36, -48, -60, -72, -84, -96, -108, -126, -144, -162, -180, -204, -234, -252, -282, -318, -354, -396, -438, -486, -534, -594, -660, -750, -798]) * -1
a_pitot, b_pitot, c_pitot = np.polyfit(calib_adc_pitot, calib_pression_pitot, 2)

### Process data
# Timestamp;Frame;Sts;Lat;Lon;Altitude;Pressure;Temperature;Acceleration X;Acceleration Y;Acceleration Z; Annex 0;Annex 1;

# time processing
tm_time = list(map(lambda x: datetime.strptime(x[0].strip("[] "), "%H:%M:%S.%f"), raw_data))
tm_time = list(map(lambda x: (x.hour * 3600 + x.minute * 60 + x.second) * 1000 + x.microsecond / 1000, tm_time))
tm_time = list(map(lambda x: x - tm_time[0], tm_time)) # msec

## Other processing
rocket_state = list(map(lambda x : (int(x[2]) & 0x7),  raw_data))

for j in range(len(rocket_state) - 1):
    if rocket_state[j] == 1 and rocket_state[j + 1] == 2:
        ascend_time = tm_time[j + 1]
tm_time = list(map(lambda x: x - ascend_time, tm_time)) # msec

gnss_valid = list(map(lambda x : (int(x[2]) & 0x8) >> 3,  raw_data))
latitude = list(map(lambda x : int(x[3]) * 1e-7,  raw_data))
longitude = list(map(lambda x : int(x[4]) * 1e-7,  raw_data))
altitude = list(map(lambda x : int(x[5])*1.0,  raw_data)) # meter
pressure = list(map(lambda x : int(x[6])/4096.0,  raw_data)) # mBar
temperature = list(map(lambda x : int(x[7])/100.0,  raw_data)) # °C
acc_x = list(map(lambda x : int(x[8])*acc_factor,  raw_data)) # g
acc_y = list(map(lambda x : int(x[9])*acc_factor,  raw_data)) # g
acc_z = list(map(lambda x : int(x[10])*acc_factor,  raw_data)) # g

## Pitot processing
pitot = list(map(lambda x : abs(int(x[11])),  raw_data))
pitot = list(map(lambda x : a_pitot*(x**2) + b_pitot*x + c_pitot,  pitot)) # mBar
# densite_air = (pression_atmospherique * masse_molaire_air) / (constante_gaz * temperature_air_kelvin)
air_density = (pressure[0]*100 * 0.02896) / (8.314 * 293.15)
# Bernoulli => speed = sqrt((2 * dynamic pressure) / density)
dynamic_pressure_pitot = list(map(lambda x,y : x-y,  pitot, pressure)) # mBar
speed_pitot = list(map(lambda x : math.sqrt((2 * x * 100) / air_density),  dynamic_pressure_pitot)) # m/s

## GNSS extrapolation
nb_gnss_point_touchdown = 20
lon_end = np.array(longitude)[gnss_valid]
lon_end = longitude[len(lon_end)-nb_gnss_point_touchdown:len(lon_end)]
lat_end = np.array(latitude)[gnss_valid]
lat_end = latitude[len(lat_end)-nb_gnss_point_touchdown:len(lat_end)]
alt_end = np.array(altitude)[gnss_valid]
alt_end = altitude[len(alt_end)-nb_gnss_point_touchdown:len(alt_end)]
time_end = np.array(tm_time)[gnss_valid]
time_end = tm_time[len(time_end)-nb_gnss_point_touchdown:len(time_end)]

# CubicSpline extrapolation 
spline_lon = CubicSpline(time_end, lon_end, extrapolate=True)
spline_lat = CubicSpline(time_end, lat_end, extrapolate=True)
spline_alt = CubicSpline(time_end, alt_end, extrapolate=True)

# Get gnss point 500 ms after last point
time_after_last_point = 500
lon_end.append(spline_lon(time_end[len(time_end)-1]+time_after_last_point))
lat_end.append(spline_lat(time_end[len(time_end)-1]+time_after_last_point))
alt_end.append(spline_alt(time_end[len(time_end)-1]+time_after_last_point))

print(alt_end[len(alt_end)-1])

### Export data
# GNSS data
kml = simplekml.Kml()
for iGnss in range(len(gnss_valid)):
    if (gnss_valid[iGnss] == 1):
        kml.newpoint(name=str(iGnss), coords=[(longitude[iGnss], latitude[iGnss], altitude[iGnss])], altitudemode = simplekml.AltitudeMode.absolute)  # lon, lat, height
kml.save("output/mse_gps_points.kml")

# Export GNSS extrapolated data
kml_extrapolated = simplekml.Kml()
for iGnss in range(len(lon_end)):
    kml_extrapolated.newpoint(name=str(iGnss), coords=[(lon_end[iGnss], lat_end[iGnss], alt_end[iGnss])], altitudemode = simplekml.AltitudeMode.absolute)  # lon, lat, height
kml_extrapolated.save("output/mse_gps_points_touchdown.kml")


### Plots
# Définition des titres et des unités pour chaque sous-graphique
titles = ['Rocket State', 'GNSS Valid', 'Latitude', 'Longitude', 'Altitude [m]', 'Pressure [mBar]',
          'Temperature [°C]', 'Acceleration X [g]', 'Acceleration Y [g]', 'Acceleration Z [g]',
          'Pitot [mBar]', 'Pitot Speed [m/s]']

units = [None, None, None, None, 'm', 'mBar', '°C', 'g', 'g', 'g', 'mBar', 'm/s']

# Définition des listes de données pour chaque sous-graphique
data_lists = [
    list(map(lambda x: const_rocket_states[x], rocket_state)),  # Rocket State
    list(map(lambda x: const_gnss_valid[x], gnss_valid)),  # GNSS Valid
    latitude,  # Latitude
    longitude,  # Longitude
    altitude,  # Altitude
    pressure,  # Pressure
    temperature,  # Temperature
    acc_x,  # Acceleration X
    acc_y,  # Acceleration Y
    acc_z,  # Acceleration Z
    pitot,  # Pitot
    speed_pitot  # Pitot Speed
]

# Tracer les courbes avec les barres
plt.figure(figsize=(12, 10))

for i in range(len(titles)):
    plt.subplot(4, 3, i + 1)
    plt.plot(tm_time, data_lists[i], label=titles[i])
    plt.xlabel('Time (ms)')
    plt.ylabel(units[i])

    # Ajouter les barres pour rocket_state passant de 1 à 2
    for j in range(len(tm_time)):
        if rocket_state[j] == 1 and rocket_state[j + 1] == 2:
            plt.axvline(x=tm_time[j + 1], color='r', linestyle='--')

    # Ajouter les barres pour rocket_state passant de 2 à 5
    for j in range(len(tm_time)):
        if rocket_state[j] == 2 and rocket_state[j + 1] == 5:
            plt.axvline(x=tm_time[j + 1], color='g', linestyle='--')

    plt.legend()

plt.tight_layout()
plt.show()

# Create new figure for GNSS data
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Display gnss points
ax.scatter(lon_end, lat_end, alt_end, c='red', marker='o', s=40, label='GNSS points')
ax.legend()
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Altitude')

plt.show()
