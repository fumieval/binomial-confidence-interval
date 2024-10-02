import matplotlib.pyplot as plt
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1])

groups = df.groupby('algorithm')

for algorithm in groups.groups.keys():
    if algorithm == 'dumb':
        continue
    if algorithm == 'wald':
        continue
    group = groups.get_group(algorithm)
    plt.plot(group['prob'], group['rate'], label=algorithm)

# Add a red line on rate = 0.95, with a label
plt.axhline(y=0.95, color='r', linestyle='--', label='95%')

plt.xlabel('Probability')
plt.ylabel('Coverage')
plt.legend()
plt.show()
