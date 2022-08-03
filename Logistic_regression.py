import pandas as pd
import statsmodels.formula.api as smf
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


# Load data
fires = pd.read_csv("fires.csv")

erc_percentiles = []
for i in fires["erc"]:
    erc_percentiles.append(stats.percentileofscore(fires["erc"], i))

fires['erc_percentiles'] = erc_percentiles
print(fires)

# Define and fit model
log_reg = smf.logit("fire ~ erc_percentiles", data=fires).fit()
# Summary of results
print(log_reg.summary())

#ercNew = pd.DataFrame({'erc': np.linspace(fires["erc"].min(), fires["erc"].max(), 1000)})
ercNew = pd.DataFrame({'erc_percentiles': np.linspace(0, 100, 1000)})
predProbs = log_reg.predict(ercNew)

plt.figure(figsize=(7, 5))
ax = plt.axes()
#ax.scatter(fires["erc"], fires["fire"], color='b', alpha=0.20)
plt.grid(linestyle='-', linewidth=0.5)
ax.scatter(ercNew, predProbs, color="black", s=1)
plt.ylim([0, 1])
ax.set_ylabel('Probability', fontsize=12)
plt.xlim([0, 100])
ax.set_xlabel('ERC-G', fontsize=12)
#plt.show()
#plt.savefig('ERC_probability.png', dpi=300)

#print(stats.percentileofscore(fires['erc'], 80))

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix

X = fires["erc"].to_numpy().reshape(-1, 1)
m2 = LogisticRegression().fit(X, fires["fire"])

print(confusion_matrix(fires["fire"], m2.predict(X)))
print("hello")
