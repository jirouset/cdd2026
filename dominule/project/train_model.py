"""
Model training workflow:

1. Call prepare_classification().           -> you get train_X, train_y, test_X, test_y
2. Train the model with train_X, train_y.   -> you get model
3. Evaluate the model with test_X, test_y.  -> you get evaluation report

Example usage:

train_X, train_y, test_X, test_y = prepare_data(df, target_col='set', feature_cols=['freq_50', 'freq_80', 'freq_100'])
trained_model = train_RF(train_X, train_y)
evaluate_model(trained_model, test_X, test_y)
"""
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# ------------------- DATA PREPARATION ---------------------------

def prepare_data(df, target_col, feature_cols):
    """
    Prepare data for binary classification using 80/20 split.

    :param df: your dataset
    :param target_col: target column name with 0s and 1s
    :param feature_cols: list of feature column names

    :return:
    train_X, train_y, test_X, test_y
    """

    X = df[feature_cols]
    y = df[target_col]

    # Split data into training and testing sets
    train_X, test_X, train_y, test_y = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    return train_X, train_y, test_X, test_y



# ------------------- MODEL TRAINING ---------------------------
"""
Choose model to be trained and tune hyperparameters!

Available models: Random Forest, KNN, logistic regression (doesn't work right now)
"""


def train_RF(X, y, n_estimators=[100, 200], max_depth=[None, 10, 20], min_samples_split=[2, 5]):
    """
    Random forest.
    It is highly robust and handles non-linear relationships well.

    :param n_estimators: number of trees in the forest
    :param max_depth: maximum depth of the tree
    :param min_samples_split: minimum number of samples required to split an internal node
    """
    rf = RandomForestClassifier(random_state=42)

    # Hyperparameters to be tuned
    param_grid = {
        'n_estimators': n_estimators,
        'max_depth': max_depth,
        'min_samples_split': min_samples_split
    }

    # Search for the best parameters
    grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X, y)

    print(f"Best Parameters: {grid_search.best_params_}")
    return grid_search.best_estimator_


def train_KNN(X, y, n_neighbors=[1, 3, 5]):
    """
    K-nearest neighbors with standardization.

    A simple, instance-based learner.
    It classifies a data point based on how its "n_neighbors" are classified.
    It is useful if your data has very distinct clusters.
    """
    knn = KNeighborsClassifier()

    # Hyperparameters to be tuned
    param_grid = {
        'n_neighbors': n_neighbors
    }

    # Search for the best parameters
    grid_search = GridSearchCV(estimator=knn, param_grid=param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X, y)

    print(f"Best Parameters: {grid_search.best_params_}")
    return grid_search.best_estimator_


def train_logistic_regression(X, y, lasso_ratio=[0.01, 0.1, 0.5, 1], c=[0.01, 0.1, 1, 10]):
    """
    Logistic regression with standardization.
    """
    pipe = Pipeline([
        ('scaler', StandardScaler()),
        # Remove penalty='elasticnet' to stop the warning
        ('lr', LogisticRegression(solver='saga', max_iter=5000, random_state=42))
    ])

    param_grid = {
        'lr__l1_ratio': lasso_ratio,
        'lr__C': c
    }

    grid_search = GridSearchCV(estimator=pipe, param_grid=param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X, y)

    print(f"Best Parameters: {grid_search.best_params_}")
    return grid_search.best_estimator_


def train_ridge_logistic(X, y, c=[0.01, 0.1, 1, 10]):
    """
    Ridge regression with standardization (both ridge and logictic).
    """

    pipe = Pipeline([
        ('scaler', StandardScaler()),
        # Leave penalty out; l1_ratio=0 defines it as Ridge (L2)
        ('lr', LogisticRegression(solver='saga', l1_ratio=0, max_iter=5000, random_state=42))
    ])

    param_grid = {
        'lr__C': c
    }

    grid_search = GridSearchCV(estimator=pipe, param_grid=param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X, y)

    print(f"Best Parameters: {grid_search.best_params_}")
    return grid_search.best_estimator_


# ------------------- MODEL EVALUATION ---------------------------

def evaluate_model(model, X_test, y_test):
    """
    Print an evaluation report including accuracy, precision, recall, and f1-score.
    """
    # Generate predictions
    y_pred = model.predict(X_test)

    # Calculate and print overall accuracy
    acc = accuracy_score(y_test, y_pred)
    print(f"Final Model Accuracy: {acc:.4f}")
    print("-" * 30)

    # Print detailed report for Precision, Recall, and F1
    print("Classification Report:")
    print(classification_report(y_test, y_pred))