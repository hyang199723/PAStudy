import requests
import json
import pandas as pd

wk_dir = '/Users/hongjianyang/Research/AirPollution'
# Get all sensors data
API_KEY = 'ED254E0B-510E-11EB-9893-42010A8001E8'
options = ['icon', 'name', 'location_type', 'latitude', 'longitude', 'altitude', 'last_seen',
          'channel_state', 'channel_flags', 'confidence', 'pm1.0', 'pm2.5',
          'pm2.5_10minute', 'pm2.5_30minute', 'pm2.5_60minute', 'pm2.5_6hour',
          'pm2.5_24hour', 'pm2.5_1week', 'pm10.0', '0.3_um_count', '0.5_um_count',
          '1.0_um_count', '2.5_um_count', '5.0_um_count', '10.0_um_count', 'humidity',
          'temperature', 'pressure', 'voc', 'ozone1', 'analog_input']
option_string = ','.join(options)
parameters = {
    "fields": option_string,
}
header = {'X-API-Key': API_KEY, 'Content-Type': 'application/json'}

response = requests.get("https://api.purpleair.com/v1/sensors", headers = header,
                        params = parameters)
res = response.json()
print(response.status_code)
#%% Process sensor readings
data = res['data']
col = res['fields']

df = pd.DataFrame(data, columns = col)
df = df[~(df['channel_state'] == 0)]
df = df[df['confidence'] >= 30]

col_to_keep = ['sensor_index', 'name', 'location_type', 'latitude', 'longitude', 'altitude',
          'channel_state', 'channel_flags', 'confidence', 'pm1.0', 'pm2.5',
          'pm2.5_10minute', 'pm2.5_30minute', 'pm2.5_60minute', 'pm2.5_6hour',
          'pm2.5_24hour', 'pm2.5_1week', 'pm10.0', '0.3_um_count', '0.5_um_count',
          '1.0_um_count', '2.5_um_count', '5.0_um_count', '10.0_um_count', 'humidity',
          'temperature', 'pressure']
df = df.loc[:, col_to_keep]

# Remove Null from longitude and latitude
df = df.dropna(subset = ['longitude', 'latitude', 'pm2.5'])

#%%
df.to_csv(wk_dir + '/Data/sensors.csv', index = False)