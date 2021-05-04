# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 22:59:25 2021

@author: lenovo
"""

import os
import sys
import time


def gpu_info():
    gpu_status = os.popen('nvidia-smi | grep %').read().split('|')
    gpu_num = (len(gpu_status) - 1) // 4
    gpu_status.pop(0)
    gpu_memory = []
    gpu_power = []
    for i in range(gpu_num):
        gpu_memory.append(int(gpu_status[4 * i + 1].split('/')[0].split('M')[0].strip()))
        gpu_power.append(int(gpu_status[4 * i].split('   ')[-1].split('/')[0].split('W')[0].strip()))

    return gpu_power, gpu_memory

def narrow_setup(interval=2):
    gpu_power, gpu_memory = gpu_info()
    i = 0
    while True:  # set waiting condition
        gpu_power, gpu_memory = gpu_info()
        gpu_power_str = ''
        gpu_memory_str = ''
        for i in range(len(gpu_power)):
            gpu_power_str += str(gpu_power[i])+' '
            gpu_memory_str += str(gpu_memory[i])+' '
            if gpu_power[i] < 20 or gpu_memory[i] < 4000:
                print('\n' + "nohup python -u train_nmt.py THEANO_FLAGS='floatX=float32,device=cuda%d,lib.cnmem=1' > nohup.out 2>&1 &"%i)
                os.system("nohup python -u train_nmt.py THEANO_FLAGS='floatX=float32,device=cuda%d,lib.cnmem=1' > nohup.out 2>&1 &"%i)
                return
        i = i % 5
        symbol = 'monitoring: ' + '>' * i + ' ' * (10 - i - 1) + '|'
        gpu_power_str = 'gpu power:%s W |' % gpu_power_str
        gpu_memory_str = 'gpu memory:%s MiB |' % gpu_memory_str
        sys.stdout.write('\r' + gpu_memory_str + ' ' + gpu_power_str + ' ' + symbol)
        sys.stdout.flush()
        time.sleep(interval)
        i += 1

if __name__ == '__main__':
    narrow_setup()