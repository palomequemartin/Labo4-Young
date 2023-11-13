import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metodos
import nidaqmx
import pyvisa as visa
import time
plt.style.use('./figuras.mplstyle')

max = 10
min = -10

def medicion_continua(duracion, fs, delay=0.1):
    
    N_puntos = int(duracion*fs)
    
    with nidaqmx.Task() as task:
        
        task.ai_channels.add_ai_voltage_chan("Dev1/ai0", max_val=max, min_val=min,
                                             terminal_config=nidaqmx.constants.TerminalConfiguration.DIFF)
        
        task.timing.cfg_samp_clk_timing(fs, sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)
        
        task.start()
        total = 0
        barra = []
        while total < N_puntos:
            time.sleep(delay)
            data_barra = task.read(number_of_samples_per_channel=nidaqmx.constants.READ_ALL_AVAILABLE)
            barra.extend(data_barra)
            total = total + len(data_barra)
        return barra



fs = 50000 #Frecuencia de muestreo
duracion = 90 #segundos

if not input(f'Enter p/ medir ({duracion}s): '):
    V_barra = medicion_continua(duracion, fs, delay=0.01)


# GUARDAR DATOS

df = pd.DataFrame()
df['Tension fotodiodo [V]'] = pd.Series(V_barra)
try:
    rm = visa.ResourceManager()
    gen_rm = rm.open_resource('USB0::0x0699::0x0346::C034165::INSTR')
    df['Frecuencia generador [Hz]'] = pd.Series(gen_rm.query_ascii_values('FREQ?'), index=[0])
    gen_rm.close()
except:
    print('Hubo error en pyvisa')
df['Duracion [s]'] = pd.Series([duracion], index=[0])
df['Frecuencia de sampleo [Hz]'] = pd.Series([fs], index=[0])

metodos.save(df, f'forzado_d{duracion}_fs{fs}', './Mediciones/Clase 3')