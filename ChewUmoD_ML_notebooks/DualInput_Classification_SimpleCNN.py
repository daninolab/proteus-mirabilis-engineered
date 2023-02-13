# -*- coding: utf-8 -*-
"""Upload Version Bacterial Colony Classifiers.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1wnn7pDVOb-0hZCtGJ_GvqqrYUZ9JdKYr

# Classifying Engineered Bacterial Patterns


This notebook is for classification of images of bacterial colonies formed by a strain engineered to respond to two different molecules by changing its colony pattern, in a different way depending on the presence of either or both molecules. With images collected at 3 concentrations of each of the 2 input molecules, there were 9 classes of images available, with ~20-30 images per class.

# Mount Drive/Imports
"""

# Commented out IPython magic to ensure Python compatibility.
import numpy as np
import sklearn
from sklearn.utils import class_weight
from sklearn.model_selection import train_test_split

import tensorflow as tf
try:
#   %tensorflow_version 2.x # enable TF 2.x in Colab
except Exception:
  pass
print('Using tensorflow version {}'.format(tf.__version__))

from tensorflow.keras import datasets, layers, models
from tensorflow.keras.preprocessing.image import load_img, ImageDataGenerator

import time
import random
import pathlib
import os
import shutil
from collections import Counter


import IPython.display as display
import matplotlib.pyplot as plt
from PIL import Image

# Mount my google drive
from google.colab import drive # will need to use verification code here
drive.mount('/content/gdrive', force_remount = True)

"""# Set up

Define some constants and functions for training.
"""

AUTOTUNE = tf.data.experimental.AUTOTUNE
DEFAULT_IMG_SIZE = 150
SHUFFLE_SIZE = 1000

def get_available_image_paths(data_root):
  #Get all image paths from my Drive
  all_image_paths = [str(path) for path in data_root.glob('**/*.tiff')]
  print('We have {} images'.format(len(all_image_paths)))
  print(all_image_paths[:2]) #get a sense of where the images are
  return all_image_paths

def get_classes_info(data_root, all_img_paths):
  #Show the classes we have ie the names of subdirectories
  label_names = sorted(item.name for item in data_root.glob('*/') if item.is_dir())

  #Determine an index for each class
  label_to_index = dict((name, index) for index, name in enumerate(label_names))
  all_labels = [label_to_index[pathlib.Path(path).parent.name] for path in all_img_paths]
  print(label_to_index)
  print(Counter(all_labels))

  #Check number of images per class
  return label_to_index, all_labels

def show(img, label):
  # im.show()
  # display.display(display.Image(image_path, width = 250, height = 250))
  plt.imshow(img)
  plt.title(label)
  plt.xticks([])
  plt.yticks([])
  print()

def get_data_generators(data_root, validation_split = 0.25, batch_size = 6, img_size=150, **kwargs):

  train_datagen = ImageDataGenerator(rescale=1./255, validation_split=validation_split, **kwargs) #NO augmentation

  # #don't have training/validation data yet
  train_dir = data_root
  validation_dir = data_root

  train_generator = train_datagen.flow_from_directory(
          train_dir,
          subset = 'training',
          target_size=(img_size, img_size),
          batch_size=batch_size,
          class_mode='sparse')

  validation_generator = train_datagen.flow_from_directory(
          validation_dir,
          subset = 'validation',
          target_size=(img_size, img_size),
          batch_size=batch_size,
          class_mode='sparse')
  
  return train_generator, validation_generator

def load_and_preprocess_image(path, img_size = DEFAULT_IMG_SIZE):
  #Can't use tensorflow preprocessing because tifs
  img = plt.imread(path)
  img = tf.image.resize(img, [img_size, img_size])
  img /= 255.0  # normalize pixels to 0,1
  return img

root_path = 'gdrive/My Drive/Danino_Lab/Patterning_Scans/cheW_umoD_Expts/Dataset_for_ML'
data_root = pathlib.Path(root_path)

"""# Preprocessing data for training

First, all available image paths are obtained, and labeled by class (the name of the subfolder they were located in). The classes are:

0 mM IPTG; 0% arabinose
0 mM IPTG; 0.1% arabinose
0 mM IPTG; 0.2% arabinose
2.5 mM IPTG; 0% arabinose
2.5 mM IPTG; 0.1% arabinose
2.5 mM IPTG; 0.2% arabinose
5 mM IPTG; 0% arabinose
5 mM IPTG; 0.1% arabinose
5 mM IPTG; 0.2% arabinose
"""

image_paths = get_available_image_paths(data_root)
label_to_index, all_labels = get_classes_info(data_root, image_paths)

"""For each image path, load the image, divide by 255, rescale to default image size (150), and store in X_full."""

X_full = []
for i, image_path in enumerate(image_paths):
  try:
    processed = load_and_preprocess_image(image_path, img_size = 150)
    X_full.append(processed)
  except Exception as e:
    print(e)
    print(image_path)
    break

y_full = all_labels

"""Convert both images & labels to tensors."""

X_full = tf.convert_to_tensor(X_full)
y_full = tf.convert_to_tensor(all_labels)

assert len(X_full) == len(y_full)
print('{} length data'.format(len(X_full)))

"""Shuffle the images, labels & split into training and validation. (Note: Did not attempt to save some for test set since there were already not many images)"""

indices = tf.range(start=0, limit=len(X_full), dtype=tf.int32)
shuffled_indices = tf.random.shuffle(indices)
X_full = tf.gather(X_full, shuffled_indices)
y_full = tf.gather(y_full, shuffled_indices)

split_index = int(0.8 * len(X_full))
X_train, y_train, X_val, y_val = X_full[:split_index], y_full[:split_index], X_full[split_index:], y_full[split_index:]

print(X_train.shape)
print(X_val.shape)
print(y_train.shape)
print(y_val.shape)

"""Compute the class weights since the classes have different numbers of images."""

class_weights = class_weight.compute_class_weight('balanced', classes = np.unique(all_labels), y = all_labels)
class_weights_dict = dict(zip(range(9), class_weights))
print(class_weights_dict)
print(label_to_index)

"""# Model Training Function

In order to train multiple models, defined a function which compiles the given model (to reset weights), augments the data if desired, creates ImageDataGenerator, and trains the model for desired number of epochs. Data augmentation used was rotation, horizontal/vertical flips, and brightness changes, which would produce images that would still be valid in real life.
"""

def train_model(model, batch_size=6, epochs=50, augmented=True, reweight_classes=True):

  model.compile(
    optimizer=tf.keras.optimizers.Adam(learning_rate=0.0001),
    loss='sparse_categorical_crossentropy',
    metrics=['accuracy']
  )
  
  datagen_args = {}
  if augmented:

    datagen_args = {
        'rotation_range': 20,
        'horizontal_flip': True,
        'vertical_flip': True,
        'brightness_range': [0.5, 1.5],
        # 'rescale': 3.0
    }
  datagen = ImageDataGenerator(**datagen_args)

  # compute quantities required for featurewise normalization
  # (std, mean, and principal components if ZCA whitening is applied)
  # datagen.fit(X_train)

  fit_generator_args = {}
  if reweight_classes:
    fit_generator_args = {
        'class_weight': class_weights_dict
    }
  # fits the model on batches with real-time data augmentation:
  history = model.fit_generator(
    datagen.flow(X_train, y_train, batch_size=batch_size),
    steps_per_epoch=len(X_train) / batch_size,
    epochs=epochs,
    validation_data = (X_val, y_val),
    **fit_generator_args 
  )

  return history

"""# CNN Model Training

Try a model with more convolutional/pooling blocks and more dense layers.
"""

# Build CNN Model
big_cnn_model = models.Sequential([
  layers.Conv2D(32, (3, 3), kernel_regularizer=tf.keras.regularizers.l2(0.01), bias_regularizer=tf.keras.regularizers.l2(0.01), input_shape=(DEFAULT_IMG_SIZE, DEFAULT_IMG_SIZE, 3)),
  layers.LeakyReLU(),
  layers.MaxPooling2D((2, 2)),
  layers.Conv2D(64, (3, 3), activation='relu'),
  layers.Conv2D(64, (3, 3), activation='relu'),
  layers.MaxPooling2D((2, 2)),
  layers.Conv2D(64, (3, 3), activation='relu'),
  layers.Conv2D(64, (3, 3), activation='relu'),
  layers.MaxPooling2D((2, 2)),
  layers.Conv2D(128, (3, 3), activation='relu'),
  layers.Conv2D(128, (3, 3), activation='relu'),
  layers.MaxPooling2D((2, 2)),
  layers.MaxPooling2D(),
  layers.Flatten(),
  layers.Dense(256, activation = 'relu'),
  layers.Dense(128, activation='relu'),
  layers.Dense(64, activation='relu'),
  layers.Dropout(rate=0.3),
  layers.Dense(9, activation='softmax')
  
])
print(big_cnn_model.summary())

histories_big_cnn_augmented = train_model(big_cnn_model, batch_size=8, epochs=50, augmented=True, reweight_classes=True)

# create a plotting function--reusing the code from HW1

def plot(history):
  
  # The history object contains results on the training and test
  # sets for each epoch
  acc = history.history['accuracy']
  val_acc = history.history['val_accuracy']
  loss = history.history['loss']
  val_loss = history.history['val_loss']

  # Get the number of epochs
  epochs = range(len(acc))

  plt.title('Training and validation accuracy')
  plt.plot(epochs, acc, color='blue', label='Train')
  plt.plot(epochs, val_acc, color='orange', label='Val')
  plt.xlabel('Epoch')
  plt.ylabel('Accuracy')
  plt.legend()

  _ = plt.figure()
  plt.title('Training and validation loss')
  plt.plot(epochs, loss, color='blue', label='Train')
  plt.plot(epochs, val_loss, color='orange', label='Val')
  plt.xlabel('Epoch')
  plt.ylabel('Loss')
  plt.legend()

plot(histories_big_cnn_augmented)

from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Tahoma']

import matplotlib as mpl
# mpl.rc('font', **{'family' : 'sans-serif', 'sans-serif' : ['Myriad Pro']})
# mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'

import seaborn as sns
sns.set_style('white')
sns.set_theme()

acc = histories_big_cnn_augmented.history['accuracy']
val_acc = histories_big_cnn_augmented.history['val_accuracy']
loss = histories_big_cnn_augmented.history['loss']
val_loss = histories_big_cnn_augmented.history['val_loss']

# Get the number of epochs
epochs = range(len(acc))
fig1 = plt.figure()
plt.title('Training and validation accuracy')
plt.plot(epochs, acc, color='navy', label='Train')
plt.plot(epochs, val_acc, color='goldenrod', label='Val')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.ylim(0, 1)
plt.legend(facecolor='white',edgecolor='gray')


fig2 = plt.figure()
plt.title('Training and validation loss')
plt.plot(epochs, loss, color='navy', label='Train')
plt.plot(epochs, val_loss, color='goldenrod', label='Val')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend(facecolor='white',edgecolor='gray')

from datetime import date
today = date.today()
date_string = today.strftime("%m%d%y")

extent = fig1.get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
save_dir = f'gdrive/My Drive/Danino_Lab/Patterning_Scans/cheW_umoD_Expts/models_and_histories/plots'
curves_path = save_dir + f'/{date_string}_simplecnn_training_acc_plots.svg'
fig1.savefig(curves_path, format='svg', bbox_inches="tight")

curves_path = save_dir + f'/{date_string}_simplecnn_training_loss_plots.svg'
fig2.savefig(curves_path, format='svg', bbox_inches="tight")
