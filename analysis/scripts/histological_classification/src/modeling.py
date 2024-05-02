import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

import tempfile

# Data manipulation
import numpy as np
import re

# Deep learning
import tensorflow as tf
from tensorflow import keras

from keras import layers
from keras.models import Model
from keras.applications import Xception, VGG16, VGG19
from tensorflow.keras.layers import GlobalAveragePooling2D, Dense, Dropout
from keras import backend as K

import seaborn as sns
import matplotlib.pyplot as plt

import sys

if tf.config.list_physical_devices('GPU'):
    print(">> Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
    strategy = tf.distribute.MirroredStrategy()
    print('>> Number of devices withing strategy: {}'.format(strategy.num_replicas_in_sync))
else:
    print(">> No GPUs found")
    strategy = tf.distribute.OneDeviceStrategy('/cpu:0')
    
from sklearn.metrics import accuracy_score, f1_score, recall_score, precision_score, roc_auc_score

def calculate_class_weights(y, multiplier=1):
    from sklearn.utils import class_weight
    import numpy as np

    y = np.where(y == 'Non Smoker', 0, 1)
    
    class_weights = class_weight.compute_class_weight(
                                        class_weight = "balanced",
                                        classes = np.unique(y),
                                        y = y                                                    
                                    )
    minority_class = np.argmin(np.bincount(y))

    class_weights[minority_class] *= multiplier

    class_weights = dict(zip(np.unique(y), class_weights))

    return class_weights

def setup_model(model_type, input_shape):
    """ Setup a tensorflow model for classification."""
                
    if model_type == "vgg19":
        base_model = VGG19(weights='imagenet', include_top=False, input_shape=input_shape)
        x = base_model.output
        x = Dropout(0.5)(x)
        x = Dense(512, activation='relu')(x)
        preds = Dense(1, activation='sigmoid')(x)
        model = Model(inputs=base_model.input, outputs=preds)

    elif model_type == "vgg16_1fc":
        base_model = VGG16(weights='imagenet', include_top=False, input_shape=input_shape)
        base_model_output = base_model.layers[-1].output
        base_model = Model(inputs=base_model.inputs, outputs=base_model_output)
        x = base_model.output
        x = tf.keras.layers.Flatten()(x)
        x = Dense(256, activation='relu')(x)
        x = Dropout(0.1)(x)
        preds = Dense(1, activation='sigmoid')(x)
        model = Model(inputs=base_model.input, outputs=preds)
        print(model.summary())
        return model
        
    elif model_type == "vgg16_avg1fc":
        base_model = VGG16(weights='imagenet', include_top=False, input_shape=input_shape)
        base_model_output = base_model.layers[-4].output
        base_model = Model(inputs=base_model.inputs, outputs=base_model_output)
        base_model = Model(inputs=base_model.inputs, outputs=GlobalAveragePooling2D()(base_model.output))
        x = base_model.output
        x = tf.keras.layers.Flatten()(x)
        x = Dense(4096, activation='relu')(x)
        x = Dense(4096, activation='relu')(x)
        preds = Dense(1, activation='sigmoid')(x)
        model = Model(inputs=base_model.input, outputs=preds)
        print(model.summary())
        return model 
    
    elif model_type == "xception_tl":
        base_model = Xception(include_top=False, weights="imagenet", input_shape=input_shape)

        # Freeze the base model
        base_model.trainable = False

        # Add new layers
        model = tf.keras.models.Sequential()
        model.add(base_model)
        model.add(layers.GlobalAveragePooling2D())
        model.add(layers.Dense(1, activation='sigmoid'))
        print(model.summary())
        return model
    
    elif model_type == "xception_ft":
        base_model = Xception(include_top=False, weights="imagenet", input_shape=input_shape)

        # Add new layers
        model = tf.keras.models.Sequential()
        model.add(base_model)
        model.add(layers.GlobalAveragePooling2D())
        model.add(layers.Dropout(0.5))
        model.add(layers.Dense(1, activation='sigmoid'))
        
        for layer in model.layers:
            layer.trainable = True
        print(model.summary())
        return model
    
    elif model_type == 'xception_tlft':
        base_model = Xception(include_top=False, weights="imagenet", input_shape=input_shape)

        # Make base_model layers non-trainable
        for layer in base_model.layers:
            layer.trainable = False

        # Add new layers
        model = tf.keras.models.Sequential()
        model.add(base_model)
        model.add(layers.GlobalAveragePooling2D())

        # Add a couple of Dense layers to increase model complexity
        model.add(layers.Dense(256, activation='relu'))
        model.add(layers.Dropout(0.5))

        model.add(layers.Dense(128, activation='relu'))
        model.add(layers.Dropout(0.5))

        model.add(layers.Dense(1, activation='sigmoid'))

        # Now we will unfreeze few top layers from the base model and make them trainable
        for layer in base_model.layers[-20:]:
            if not isinstance(layer, layers.BatchNormalization):
                layer.trainable = True
        
        print(model.summary())
        return model
    
    elif model_type == 'xception_ftv2':
        base_model = Xception(include_top=False, weights="imagenet", input_shape=input_shape)

        # Add new layers
        model = tf.keras.models.Sequential()
        model.add(base_model)
        model.add(layers.GlobalAveragePooling2D())

        # Add a couple of Dense layers to increase model complexity
        model.add(layers.Dense(256, activation='relu'))
        model.add(layers.Dropout(0.2))

        model.add(layers.Dense(128, activation='relu'))
        model.add(layers.Dropout(0.2))

        model.add(layers.Dense(1, activation='sigmoid'))
        
        print(model.summary())
        return model
    else:
        raise ValueError(f"Unsupported model_type {model_type}")
    return model

def setup_batch_generator(data, image_size, batch_size, subset, augmentation=None, weights=None):
    if subset == 'train':
        if augmentation:
            print('[+] Keras augmentation activated.')
            datagen_train = keras.preprocessing.image.ImageDataGenerator(
                rescale=1./255,
                rotation_range=20,
                width_shift_range=0.1,
                height_shift_range=0.1,
                # shear_range=0.1,
                # zoom_range=0.1,
                horizontal_flip=True,
                vertical_flip=True,
                fill_mode='wrap'
            )
        else:
            datagen_train = keras.preprocessing.image.ImageDataGenerator(
                rescale=1./255,
            )

        # Create the train generator
        train_generator = datagen_train.flow_from_dataframe(
            data,
            directory=None,
            x_col="tile_path",
            y_col="smoker_status",
            weight_col=weights,
            target_size=(image_size, image_size),
            batch_size=batch_size,
            class_mode="binary"
        )
        
        return train_generator
    
    elif subset == 'valid':
        datagen_valid = keras.preprocessing.image.ImageDataGenerator(
            rescale=1./255
        )

        # Create the validation generator
        validation_generator = datagen_valid.flow_from_dataframe(
            data,
            directory=None,
            x_col="tile_path",
            y_col="smoker_status",
            target_size=(image_size, image_size),
            batch_size=batch_size,
            class_mode="binary"
        )
        
        return validation_generator
    
    elif subset == 'test':
        test_datagen = keras.preprocessing.image.ImageDataGenerator(rescale=1./255)

        test_generator = test_datagen.flow_from_dataframe(
            data,
            directory=None,
            x_col="tile_path",
            y_col="smoker_status",
            target_size=(image_size, image_size),
            batch_size=batch_size,
            class_mode="binary",
            shuffle=False  # Important to keep data in the same order as labels
        )

        return test_generator
    
def keras_f1_score(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    recall = true_positives / (possible_positives + K.epsilon())
    f1_val = 2*(precision*recall)/(precision+recall+K.epsilon())
    return f1_val

def plot_roc_curve(y_test, y_pred, title, model_name='Keras'):
    from sklearn.metrics import roc_curve, auc

    # Get the ROC curve
    fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred)
    auc_keras = auc(fpr_keras, tpr_keras)

    # Plot ROC curve
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label='{} (area = {:.3f})'.format(model_name, auc_keras))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    
    file_suffix = title + '_.png'
    with tempfile.NamedTemporaryFile(suffix=file_suffix) as temp_file:
        plt.savefig(temp_file.name, format='png', dpi=80, bbox_inches="tight")
        temp_file.flush()
        mlflow.log_artifact(temp_file.name, 'figure__ROC')
    plt.close('all')

    return
