# -*- coding: utf-8 -*-
"""3_Visualize_Training_Results.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1EV9v88Bb2bkTED9uXCsXy_ybNDyKGScV
"""

import json
import matplotlib.pyplot as plt

results = None
with open('results.json', 'r') as jsonfile:
    results = json.load(jsonfile)
results

# one of 'accuracy', 'top_k_categorical_accuracy', 'val_accuracy', 'val_top_k_categorical_accuracy'
training_key ='val_top_k_categorical_accuracy'

fig, ax = plt.subplots(figsize=(7, 5))
for model_key, data in results.items():
    if model_key.endswith('pretrained'):
        plt.plot(
            range(1, len(data['training_history'][training_key]) + 1),
            data['training_history'][training_key],
            label = model_key.split('_')[0].title()
        )
for x in [10, 20, 30, 40, 50]:
    plt.axvline(x=x, color="#d3d3d3", linestyle='--')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.title('Top-3 Validation accuracy for deep CNN classifiers')
plt.legend()
plt.show()

# ax.set_rasterized(True)
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig('val-top-k-accuracy-plot.eps', format='eps', bbox_inches=extent.expanded(2.0, 2.0))

