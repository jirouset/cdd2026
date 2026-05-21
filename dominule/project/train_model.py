"""
Model training workflow fo BALANCED datasets:

1. Call prepare_classification().           -> you get train_X, train_y, test_X, test_y
2. Train the model with train_X, train_y.   -> you get model
3. Evaluate the model with test_X, test_y.  -> you get evaluation report

Example usage:

train_X, train_y, test_X, test_y = prepare_data(df, target_col='set', feature_cols=['freq_50', 'freq_80', 'freq_100'])
trained_model = train_RF(train_X, train_y)
evaluate_model(trained_model, test_X, test_y)
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, roc_curve, auc
from sklearn.model_selection import train_test_split, GridSearchCV, permutation_test_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier


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

Available classification models: Random Forest, KNN, logistic regression, ridge regression (classification).

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
        ('lr', LogisticRegression(solver='saga', max_iter=5000, random_state=42))
    ])

    param_grid = {
        'lr__l1_ratio': lasso_ratio, # elastic net
        'lr__C': c
    }

    grid_search = GridSearchCV(estimator=pipe, param_grid=param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X, y)

    print(f"Best Parameters: {grid_search.best_params_}")
    return grid_search.best_estimator_


def train_ridge_regression(X, y, c=[0.01, 0.1, 1, 10]):
    """
    Ridge regression with standardization.
    Note: Model for classification task (uses elastic net penalty in LogisticRegression).
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


def train_XGB(X, y, n_estimators=[100, 200], max_depth=[3, 6], learning_rate=[0.05, 0.1, 0.3]):
    """
    XGBoost gradient boosting classifier.
    Handles non-linear relationships and feature interactions well;
    often outperforms Random Forest on tabular data.

    :param n_estimators: number of boosting rounds
    :param max_depth: maximum tree depth
    :param learning_rate: step size shrinkage to prevent overfitting
    """
    xgb = XGBClassifier(
        random_state=42,
        eval_metric="logloss",
        verbosity=0,
    )

    param_grid = {
        "n_estimators": n_estimators,
        "max_depth": max_depth,
        "learning_rate": learning_rate,
    }

    grid_search = GridSearchCV(estimator=xgb, param_grid=param_grid, cv=5, scoring="accuracy")
    grid_search.fit(X, y)

    print(f"Best Parameters: {grid_search.best_params_}")
    return grid_search.best_estimator_


# ------------------- MODEL EVALUATION ---------------------------

def evaluate_model(model, X_test, y_test):
    """
    Print an evaluation report including accuracy, precision, recall, and f1-score.
    """
    y_pred = model.predict(X_test)

    acc = accuracy_score(y_test, y_pred)
    print(f"Final Model Accuracy: {acc:.4f}")
    print("-" * 30)

    print("Classification Report:")
    print(classification_report(y_test, y_pred))


def plot_roc_auc(models: dict, X_test, y_test):
    """
    Plot ROC curves for multiple fitted models with AUC scores in the legend.

    :param models: dict of {label: fitted_model}, e.g. {"RF": rf_model, "XGB": xgb_model}
    :param X_test: test features
    :param y_test: true binary labels
    """
    fig, ax = plt.subplots(figsize=(7, 6))

    for name, model in models.items():
        y_score = model.predict_proba(X_test)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, y_score)
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, lw=2, label=f"{name}  (AUC = {roc_auc:.3f})")

    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", lw=1)
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curves")
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.show()


def plot_feature_importance(model, feature_names, top_k: int = 10):
    """
    Plot the top k most important features for a fitted model.

    Supports RandomForest and XGBoost (feature_importances_) and
    Pipeline-wrapped LogisticRegression / Ridge (uses |coef_|).

    :param model: fitted model or Pipeline
    :param feature_names: list of feature names matching the training columns
    :param top_k: number of top features to show (default 10)
    """
    # Extract importances from model or final Pipeline step
    if isinstance(model, Pipeline):
        estimator = model[-1]
        importances = np.abs(estimator.coef_[0])
        importance_label = "Absolute coefficient"
    else:
        importances = model.feature_importances_
        importance_label = "Feature importance"

    feature_names = list(feature_names)
    indices = np.argsort(importances)[::-1][:top_k]
    top_names = [str(feature_names[i]) for i in indices]
    top_vals = importances[indices]

    fig, ax = plt.subplots(figsize=(8, max(4, top_k * 0.4)))
    ax.barh(range(top_k), top_vals[::-1], color="#5b8db8")
    ax.set_yticks(range(top_k))
    ax.set_yticklabels(top_names[::-1])
    ax.set_xlabel(importance_label)
    ax.set_title(f"Top {top_k} Feature Importances")
    plt.tight_layout()
    plt.show()


def permutation_test(
    model,
    X,
    y,
    n_permutations: int = 500,
    cv: int = 5,
    scoring: str = "accuracy",
    model_name: str = None,
    random_state: int = 42,
) -> tuple:
    """
    Permutation test: checks whether the model performs significantly better than chance.

    The model is re-trained (cloned) n_permutations times on randomly shuffled labels.
    The p-value is the fraction of permuted scores that are >= the true CV score.

    Prints a text report and shows a histogram of permuted scores with the true score marked.

    :param model:          Fitted model returned by train_RF / train_KNN / etc.
                           It is cloned internally, so the original object is untouched.
    :param X:              Feature matrix (use train_X).
    :param y:              True labels (use train_y).
    :param n_permutations: Number of label permutations (default 500; use 1000 for publication).
    :param cv:             Number of cross-validation folds used for scoring (default 5).
    :param scoring:        Sklearn scoring metric (default 'accuracy').
    :param model_name:     Display name for the report title. Inferred from class name if None.
    :param random_state:   Random seed for reproducibility.
    :return:               (true_score, permuted_scores, p_value)
    """
    if model_name is None:
        # Unwrap Pipeline to get the final estimator's name
        est = model[-1] if isinstance(model, Pipeline) else model
        model_name = type(est).__name__

    true_score, perm_scores, p_value = permutation_test_score(
        model, X, y,
        n_permutations=n_permutations,
        cv=cv,
        scoring=scoring,
        random_state=random_state,
        n_jobs=-1,
    )

    # Significance label
    if p_value < 0.001:
        sig_label = "*** (p < 0.001)"
    elif p_value < 0.01:
        sig_label = "** (p < 0.01)"
    elif p_value < 0.05:
        sig_label = "* (p < 0.05)"
    else:
        sig_label = "n.s. (p ≥ 0.05)"

    # Plot
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(perm_scores, bins=30, color="#5b8db8", edgecolor="white", alpha=0.85,
            label=f"Permuted scores (n={n_permutations})")
    ax.axvline(true_score, color="#e74c3c", linewidth=2.5,
               label=f"True CV {scoring} = {true_score:.3f}")
    ax.axvline(np.mean(perm_scores), color="grey", linewidth=1.2, linestyle="--",
               label=f"Mean permuted = {np.mean(perm_scores):.3f}")
    ax.set_xlabel(scoring.capitalize())
    ax.set_ylabel("Count")
    ax.set_title(f"Permutation Test — {model_name}\np = {p_value:.4f}  {sig_label}")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.show()

    # Text report
    interpretation = (
        "Model IS significantly better than chance."
        if p_value < 0.05
        else "Model is NOT significantly better than chance."
    )
    print(f"Permutation Test — {model_name}")
    print("-" * 45)
    print(f"  True CV {scoring:<12}: {true_score:.4f}")
    print(f"  Mean permuted     : {np.mean(perm_scores):.4f}")
    print(f"  Std  permuted     : {np.std(perm_scores):.4f}")
    print(f"  p-value           : {p_value:.4f}  {sig_label}")
    print(f"  → {interpretation}")

    return true_score, perm_scores, p_value