# difference.py
# pip install difference
from differences import ATTgt, simulate_data
from memory_profiler import profile
import time

# almost 2500 rows
df = simulate_data(nentity=313)

start_time = time.time()

@profile
def diffe_rence():
	att_gt = ATTgt(data=df, cohort_name='cohort')
	att_gt.fit(formula='y')

if __name__ == '__main__':
	diffe_rence()
	end_time = time.time()
	count = end_time - start_time
	print(f'Execution time (difference): {count}')