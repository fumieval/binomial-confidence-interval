import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('result.csv')

groups = df.groupby('algorithm')

for algorithm in groups.groups.keys():
    if algorithm == 'dumb':
        continue
    if algorithm == 'wald':
        continue
    group = groups.get_group(algorithm)
    plt.plot(group['prob'], group['avgWidth'], label=algorithm)

plt.xlabel('Probability')
plt.ylabel('Average Width')
plt.legend()
plt.show()
