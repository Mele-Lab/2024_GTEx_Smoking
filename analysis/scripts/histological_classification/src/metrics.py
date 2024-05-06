import seaborn as sns
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, f1_score, recall_score, precision_score, roc_auc_score, roc_curve, auc
import pandas as pd

def plot_cm(df, col_true='true', col_pred='pred', title=None, path=None):
    # Generate the confusion matrix
    cm = confusion_matrix(df[col_true], df[col_pred])

    # Create a DataFrame from the confusion matrix
    cm_df = pd.DataFrame(cm, index=['non-smoker', 'smoker'], columns=['non-smoker', 'smoker'])

    # Create a heatmap from the DataFrame
    sns.heatmap(cm_df, annot=True, cmap='Blues', fmt='g')
    plt.title('Confusion Matrix')
    plt.ylabel('Actual')
    plt.xlabel('Predicted')

    # Plot/save figure
    if title:
        plt.title(title)

    plt.savefig("confusion_matrix.png", format='png', dpi=80, bbox_inches="tight")
    plt.close('all')

def plot_roc_curve(y_test, y_pred, title, model_name='Keras'):
    # Get the ROC curve
    fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred)
    auc_keras = auc(fpr_keras, tpr_keras)

    # Plot/save ROC curve
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label='{} (area = {:.3f})'.format(model_name, auc_keras))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    
    plt.savefig("roc_curve.png", format='png', dpi=80, bbox_inches="tight")
    plt.close('all')

    return

def classification_metrics(y_true, y_pred, y_pred_prob, name):
    metrics = {
        f'{name}_accuracy': accuracy_score(y_true, y_pred),
        f'{name}_f1_score': f1_score(y_true, y_pred),
        f'{name}_f1_score_micro': f1_score(y_true, y_pred, average='micro'),
        f'{name}_f1_score_macro': f1_score(y_true, y_pred, average='macro'),
        f'{name}_recall': recall_score(y_true, y_pred),
        f'{name}_precision': precision_score(y_true, y_pred),
        f'{name}_auc': roc_auc_score(y_true, y_pred_prob),
    }
    
    return metrics