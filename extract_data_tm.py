import csv
import matplotlib.pyplot as plt
import numpy as np
import simplekml
import math
from time import time
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D
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
calib_pression_pitot = np.array([1.00, 1.04, 1.07, 1.11, 1.16, 1.20, 1.25, 1.31, 1.37, 1.43, 1.50, 1.58, 1.67, 1.77, 1.88, 2.00, 2.15, 2.31, 2.50, 2.73, 3.00, 3.34, 3.76, 4.29, 5.01, 6.01, 7.51]) * 1000
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

# Find time index corresponding to rocket state
st_ascend_index = 0
st_descend_index = 0
for j in range(len(tm_time)):
    if rocket_state[j] == 1 and rocket_state[j + 1] == 2:
        st_ascend_index = tm_time[j + 1]
    if rocket_state[j] == 2 and rocket_state[j + 1] == 5:
        st_descend_index = tm_time[j + 1]

gnss_valid = list(map(lambda x : (int(x[2]) & 0x8) >> 3,  raw_data))
latitude = list(map(lambda x : int(x[3]) * 1e-7,  raw_data))
longitude = list(map(lambda x : int(x[4]) * 1e-7,  raw_data))
altitude = list(map(lambda x : int(x[5])*1.0,  raw_data)) # meter
pressure = list(map(lambda x : int(x[6])/4096.0,  raw_data)) # mBar
temperature = list(map(lambda x : int(x[7])/100.0,  raw_data)) # °C
acc_x = list(map(lambda x : int(x[8])*acc_factor,  raw_data)) # g
acc_y = list(map(lambda x : int(x[9])*acc_factor,  raw_data)) # g
acc_z = list(map(lambda x : int(x[10])*acc_factor,  raw_data)) # g

## Barometer altitude
# altitude_baro = list(map(lambda x,y: (288.15 / 0.0065) * (1 - (x*100 / 101325) ** ( (9.81 * 0.0065) / 8.314)), pressure, temperature))
altitude_baro = list(map(lambda x,y: -(math.log(x/pressure[0]) * 8.31432 * (y+273.15))/(9.81*0.0289644), pressure, temperature))
max_alt_ind = altitude_baro.index(max(altitude_baro))
# altitude_baro = list(map(lambda x,y: -(math.log(x/pressure[0]) * 8.31432 * (y+273.15))/(9.81*0.0289644), pressure, temperature))

speed_baro = [0]
for i in range(0, len(tm_time)-1, 1):
    sp = abs((altitude_baro[i+1] - altitude_baro[i]) / (tm_time[i+1] - tm_time[i])) * 1000
    speed_baro.append(sp)

max_speed_ind = speed_baro.index(max(speed_baro))

## Pitot processing
pitot = list(map(lambda x : abs(int(x[11])),  raw_data))
pitot = list(map(lambda x : a_pitot*(x**2) + b_pitot*x + c_pitot,  pitot)) # mBar
# densite_air = (pression_atmospherique * masse_molaire_air) / (constante_gaz * temperature_air_kelvin)
air_density = (pressure[0]*100 * 0.02896) / (8.314 * 293.15)
# Bernoulli => speed = sqrt((2 * dynamic pressure) / density)
static_pressure_diff = pitot[0] - pressure[0]
dynamic_pressure_pitot = list(map(lambda x,y : x-y-static_pressure_diff if (x-y-static_pressure_diff >= 0) else 0,  pitot, pressure)) # mBar
speed_pitot = list(map(lambda x : math.sqrt((2 * x * 100) / air_density),  dynamic_pressure_pitot)) # m/s

## GNSS extrapolation
nb_gnss_point_touchdown = 20
# lon_end = np.array(longitude)[gnss_valid]
lon_end = longitude[len(longitude)-nb_gnss_point_touchdown:len(longitude)]
# lat_end = np.array(latitude)[gnss_valid]
lat_end = latitude[len(latitude)-nb_gnss_point_touchdown:len(latitude)]
# alt_end = np.array(altitude)[gnss_valid]
alt_end = altitude[len(altitude)-nb_gnss_point_touchdown:len(altitude)]
# time_end = np.array(tm_time)[gnss_valid]
time_end = tm_time[len(tm_time)-nb_gnss_point_touchdown:len(tm_time)]

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
titles = ['Rocket State', 'GNSS Valid', 'Latitude', 'Longitude', 'Altitude [m]', 'Pressure [mBar]',
          'Temperature [°C]', 'Acceleration X [g]', 'Acceleration Y [g]', 'Acceleration Z [g]',
          'Pitot [mBar]', 'Pitot Speed [m/s]', 'Altitude barometer']
units = [None, None, None, None, 'm', 'mBar', '°C', 'g', 'g', 'g', 'mBar', 'm/s', 'm']
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
    speed_pitot,  # Pitot Speed
    altitude_baro
]

def plot_rocket_states(ax):
    ax.axvline(x=st_ascend_index, color='r', linestyle='--')
    ax.axvline(x=st_descend_index, color='g', linestyle='--')
    ax.axvline(x=tm_time[max_alt_ind], color='g', linestyle=':')

fig_0 = plt.figure(figsize=(12, 10))

for i in range(len(titles)):
    plt.subplot(4, 4, i + 1)
    plt.plot(tm_time, data_lists[i], label=titles[i])
    plt.xlabel('Time (ms)')
    plt.ylabel(units[i])
    plot_rocket_states(plt)
    plt.legend()

plt.suptitle('MSE telemetry', fontsize=16, fontweight='bold', color='blue')
plt.tight_layout()

## Combined plots
fig_1 = plt.figure(figsize=(12, 12))

ax1 = plt.subplot(2, 1, 1)
ax1.plot(tm_time, list(map(lambda x, y: x-altitude[0] if y == 1 else 0, altitude, gnss_valid)), label='GNSS "valid" altitude (relative to surface)')
ax1.plot(tm_time, altitude_baro, label='Barometer altitude (relative to surface)')
ax1.scatter(tm_time[max_alt_ind], altitude_baro[max_alt_ind], marker='+', s=200)
ax1.text(tm_time[max_alt_ind], altitude_baro[max_alt_ind], f'Max alt: {altitude_baro[max_alt_ind]:.0f} m @ {tm_time[max_alt_ind]/1000.0:.2f} s', ha='left', va='bottom', color='blue')
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('m')
plot_rocket_states(ax1)
plt.legend()

ax1 = plt.subplot(2, 1, 2)
ax1.plot(tm_time, speed_baro, label='Barometer speed')
ax1.plot(tm_time, speed_pitot, label='Pitot speed')
ax1.scatter(tm_time[max_speed_ind], speed_baro[max_speed_ind], marker='+', s=200)
ax1.text(tm_time[max_speed_ind], speed_baro[max_speed_ind], f'Max speed: {speed_baro[max_speed_ind]:.2f} m/s @ {tm_time[max_speed_ind]/1000.0:.2f} s', ha='left', va='bottom', color='blue')
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('m/s')
plot_rocket_states(ax1)
plt.legend()

plt.tight_layout()

### GNSS plots
fig_2 = plt.figure(figsize=(10, 8))
ax = fig_2.add_subplot(111, projection='3d')

ax.scatter(lon_end, lat_end, alt_end, c='red', marker='o', s=40, label='GNSS Points')
ax.legend()
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Altitude')

plt.suptitle('GNSS data', fontsize=16, fontweight='bold', color='blue')
plt.tight_layout()
plt.show()
