# Genomic Sequence Classification for Transcription Factor Binding Site Prediction

This Python Jupyter notebook implements a Convolutional Neural Network (CNN) using PyTorch to predict transcription factor (TF) binding sites in DNA sequences by identifying motifs. It includes data loading, preprocessing, model training, evaluation, and visualization. Below is a summary of the key steps:

**Library Imports**  
Imports essential Python libraries: PyTorch (for CNN implementation), NumPy (for array operations), scikit-learn (for data splitting and evaluation), Matplotlib (for loss plots), and Seaborn (for confusion matrix visualization).

**Data Loading**  
Reads DNA sequences from `sequences.txt` (2000 sequences, each 50 base pairs) and binary labels (0 or 1) from `labels.txt`, where 1 indicates a sequence with a potential TF binding motif and 0 indicates a non-binding sequence.

**Data Preprocessing**  
- Applies one-hot encoding to DNA sequences (A, C, G, T → 4D vectors), resulting in a tensor of shape `[num_sequences, 4, 50]`.  
- Splits data into training (60%, 1200 sequences), validation (20%, 400 sequences), and test (20%, 400 sequences) sets using scikit-learn’s `train_test_split`.

**CNN Model Definition**  
Defines a `GenomicCNN` class with:  
- Two 1D convolutional layers (4 → 16 channels, 16 → 32 channels, kernel_size=3, padding=1).  
- ReLU activations and max-pooling layers (kernel_size=2) to capture motif patterns.  
- A fully connected layer and sigmoid activation for binary classification (TF binding vs. non-binding).

**Model Training**  
- Trains the CNN using Adam optimizer (learning rate: 0.001) and Binary Cross-Entropy loss.  
- Uses mini-batch gradient descent (batch size: 32) for up to 1000 epochs with early stopping (patience: 10 epochs).  
- Saves training and validation loss plots every 10 epochs as `loss_epoch_<epoch>.png`.  
- Saves the trained model weights to `genomic_cnn_model.pth`.

**Model Evaluation**  
- Evaluates the model on the test set, computing accuracy (99.25%) and a classification report (precision, recall, F1-score).  
- Generates a confusion matrix heatmap to visualize true vs. predicted labels, saved as `confusion_matrix.png`.

**Key Outputs**  
- **Model File**: `genomic_cnn_model.pth` (trained model weights).  
- **Plots**: Loss curves (`loss_epoch_<epoch>.png`) and confusion matrix (`confusion_matrix.png`).  
- **Metrics**: Accuracy and classification report printed to the console.

This notebook provides a robust pipeline for predicting transcription factor binding sites in DNA sequences, enabling motif discovery and supporting computational genomics research.