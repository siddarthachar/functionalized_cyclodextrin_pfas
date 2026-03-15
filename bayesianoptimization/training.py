import torch
import gpytorch
from gpytorch.mlls import ExactMarginalLogLikelihood

def train_model(model, mll, train_x, train_y, max_iter=1000, lr=0.1, print_freq=500):
    """
    Train a GPyTorch Gaussian Process model using ExactMarginalLogLikelihood.

    Args:
        model: GPyTorch ExactGP model
        mll: Marginal log likelihood (usually ExactMarginalLogLikelihood)
        train_x (torch.Tensor): Input features
        train_y (torch.Tensor): Target values
        max_iter (int): Number of training iterations
        lr (float): Learning rate
        print_freq (int): Frequency of loss printing

    Returns:
        model: Trained model
    """
    model.train()
    mll.likelihood.train()

    if train_y.dim() > 1:
        train_y = train_y.squeeze()

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    # optimizer = torch.optim.LBFGS(model.parameters(), lr=lr, line_search_fn='strong_wolfe')

    for i in range(max_iter):
        
        optimizer.zero_grad()
        output = model(train_x)
        loss = -mll(output, train_y)
        loss.backward()
        if (i + 1) % print_freq == 0 or i == 0:
            print(f"[Iter {i+1}/{max_iter}] Loss: {loss.item():.4f}")
        optimizer.step()
    return model
