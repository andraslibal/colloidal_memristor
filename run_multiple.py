import os
import time
import psutil

relax_durations = [1000000, 2000000, 3000000, 4000000, 5000000]
temperatures = [2]

for rd in relax_durations:
	print('relax duration = {}'.format(rd))

	for t in temperatures:
		print(f"t = {t}")

		prefix = f'relax_d{rd}'

		f = open("parameters_base.txt", "rt")
		data = f.read()
		data = data.replace('basename', f'basename {prefix}')
		data = data.replace('relax_duration', f'relax_duration {rd}')
		data = data.replace('temperature', f'temperature {t}')
		f.close()

		f = open("parameters.txt", "wt")
		f.write(data)
		f.close()

		os.system(f"./memristor parameters.txt > ./particles/log_{prefix} &")

		time.sleep(5)

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