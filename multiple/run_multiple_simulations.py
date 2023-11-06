import os
import time
import psutil
import sys

def wait():
    elapsed = 0

    while (1):
        tmp = psutil.cpu_percent()

        if (tmp > 95.0):
            print('\r --- waiting for a free CPU ... elapsed time: {} sec ... cpu percent: {}%'.format(elapsed*5, tmp), end='')
            elapsed += 1
            time.sleep(5)
        else:
            if (elapsed):
                print()
            break

for temp in [0.5, 1, 2, 4, 8, 16]:
    print(f'temperature = {temp}')

    prefix = f't{temp}'

    name = f"{prefix}.mvi"

    f = open("parameters_base.txt", "rt")
    data = f.read()
    data = data.replace('basename', f'basename {prefix}')
    data = data.replace('temperature', f'temperature {temp}')
    f.close()

    f = open("parameters.txt", "wt")
    f.write(data)
    f.close()

    os.system(f"./memristor ./parameters.txt > ./particles/log_{prefix} &")
    time.sleep(5)
    wait()
