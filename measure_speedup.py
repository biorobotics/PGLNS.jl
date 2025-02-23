import numpy as np
import os

def get_cost_vs_time(tour_file):
  with open(tour_file, 'r') as f:
    lines = f.readlines()

  instance_name = lines[0][lines[0].find(":") + 2:-1]

  cost_vs_time = []
  started_tour_history = False
  for line_idx, line in enumerate(lines):
    if 'Tour History' in line:
      started_tour_history = True
      continue
    if not started_tour_history:
      continue

    timestamp = float(line[line.find('(') + 1:line.find(',')])
    cost = int(line[line.find(']') + 3:line.find(')')])
    cost_vs_time.append((timestamp, cost))

  return instance_name, np.array(cost_vs_time)

folders = ["solutions_1threads/", "solutions_8threads/"]
cost_vs_time_arrs = []
for folder in folders:
  cost_vs_time_arrs.append(dict())

  for filename in os.listdir(folder):
    if 'baf' not in filename:
      continue
    instance_name, cost_vs_time = get_cost_vs_time(folder + "/" + filename)
    cost_vs_time_arrs[-1][instance_name] = cost_vs_time

speedups = []
percent_differences_in_cost = []
num_instances = 0
num_instances_with_speedup = 0
for k in cost_vs_time_arrs[0]:
  cost_vs_time_1thread = cost_vs_time_arrs[0][k]
  if k not in cost_vs_time_arrs[1]:
    continue
  cost_vs_time_8threads = cost_vs_time_arrs[1][k]

  final_cost_1thread = cost_vs_time_1thread[-1, 1]
  final_time_1thread = cost_vs_time_1thread[-1, 0]

  final_cost_8threads = cost_vs_time_8threads[-1, 1]
  final_time_8threads = cost_vs_time_8threads[-1, 0]

  reach_idx_8threads = np.where(cost_vs_time_8threads[:, 1] <= final_cost_1thread)[0]
  if len(reach_idx_8threads):
    reach_time_8threads = cost_vs_time_8threads[reach_idx_8threads[0], 0]
    speedups.append(final_time_1thread/reach_time_8threads)
    num_instances_with_speedup += 1
  else:
    print(final_cost_8threads, final_cost_1thread)
  num_instances += 1

  percent_differences_in_cost.append((final_cost_8threads - final_cost_1thread)/final_cost_1thread)

print('Min, median, max speedup')
print(np.amin(speedups), np.median(speedups), np.amax(speedups))
print('Number of instances with speedup, number of instances, and percentage of instances with speedup')
print(num_instances_with_speedup, num_instances, num_instances_with_speedup/num_instances*100)

print('Min, median, max percent difference of 8thread cost from 1thread cost')
print(np.amin(percent_differences_in_cost), np.median(percent_differences_in_cost), np.amax(percent_differences_in_cost))

costs_1thread = [v[-1, 1] for v in cost_vs_time_arrs[0].values()]
costs_8threads = [v[-1, 1] for v in cost_vs_time_arrs[1].values()]

print('Min, median, max cost (1 thread)')
print(np.amin(costs_1thread), np.median(costs_1thread), np.amax(costs_1thread))

print('Min, median, max cost (8 threads)')
print(np.amin(costs_8threads), np.median(costs_8threads), np.amax(costs_8threads))
