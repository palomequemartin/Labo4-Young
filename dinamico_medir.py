import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metodos
import nidaqmx
import time
plt.style.use('./figuras.mplstyle')


def medicion_continua(duracion, fs, delay=0.1):
    
    N_puntos = int(duracion*fs)
    
    with nidaqmx.Task() as task:
        
        task.ai_channels.add_ai_voltage_chan("Dev1/ai0", max_val=10, min_val=-10,
                                             terminal_config=nidaqmx.constants.TerminalConfiguration.DIFF)
        
        task.timing.cfg_samp_clk_timing(fs, sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)
        
        task.start()
        total = 0
        data = []
        tiempo = []
        datos_hasta_tiempo = []
        t0 = time.time()
        while total < N_puntos:
            time.sleep(delay)
            datos = task.read(number_of_samples_per_channel=nidaqmx.constants.READ_ALL_AVAILABLE)
            data.extend(datos)
            total = total + len(datos)
            tiempo.append(time.time() - t0)
            datos_hasta_tiempo.append(len(datos))
        return data, tiempo, datos_hasta_tiempo



fs = 50000 #Frecuencia de muestreo
duracion = 420 #segundos

if not input(f'Enter p/ medir ({duracion}s): '):
    V, t, datos_hasta_tiempo = medicion_continua(duracion, fs)


# GUARDAR DATOS

df = pd.DataFrame()
df['Tension fotodiodo [V]'] = pd.Series(V)
df['Tiempo [s]'] = pd.Series(t)
df['Datos hasta tiempo'] = pd.Series(datos_hasta_tiempo)
df['Duracion [s]'] = pd.Series([duracion], index=[0])
df['Frecuencia de sampleo [Hz]'] = pd.Series([fs], index=[0])

metodos.save(df, f'dinamico_d{duracion}_fs{fs}', './Mediciones/Clase 3')