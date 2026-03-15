import torch
import gpytorch
from botorch.models.gpytorch import GPyTorchModel

class AdditiveGPModel_botorch(gpytorch.models.ExactGP, GPyTorchModel):
    def __init__(self, train_x, train_y, likelihood,
                 primary_dim, secondary_dim,
                 nu=2.5,
                 init_lengthscale_primary=5,
                 init_lengthscale_secondary=5):
        
        super(AdditiveGPModel_botorch, self).__init__(train_x, train_y, likelihood)
        # self.num_outputs = 1 
        self.mean_module = gpytorch.means.ConstantMean()
        self.primary_dim = primary_dim
        self.secondary_dim = secondary_dim
        self.covar_module_primary = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=nu, ard_num_dims=primary_dim)
        )
        self.covar_module_primary.base_kernel.lengthscale = init_lengthscale_primary
        self.covar_module_secondary = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=nu, ard_num_dims=secondary_dim)
        )
        self.covar_module_secondary.base_kernel.lengthscale = init_lengthscale_secondary
        # self.covar_module_primary.outputscale = 1.0  # or some appropriate value
        # self.covar_module_secondary.outputscale = 1.0
    @property
    def num_outputs(self):
        return 1

    def forward(self, x):
        mean_x = self.mean_module(x)
        x_primary = x[..., :self.primary_dim]
        x_secondary = x[..., self.primary_dim:self.primary_dim + self.secondary_dim]

        # x_primary = x[:, :self.primary_dim]
        # x_secondary = x[:, self.primary_dim:self.primary_dim + self.secondary_dim]
        covar_primary = self.covar_module_primary(x_primary)
        covar_secondary = self.covar_module_secondary(x_secondary)
        covar_total = covar_primary + covar_secondary
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_total)


# defining kernel for GP - check below for other models
class AdditiveGPModel_gpytorch(gpytorch.models.ExactGP): # this is just for gpytorch models. The one above is for botorch
    def __init__(self, train_x, train_y, likelihood,
                 primary_dim, secondary_dim,
                 nu=2.5,
                 init_lengthscale_primary=5,
                 init_lengthscale_secondary=5):
        
        super(AdditiveGPModel_gpytorch, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()

        self.primary_dim = primary_dim
        self.secondary_dim = secondary_dim

        self.covar_module_primary = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=nu, ard_num_dims=primary_dim)
        )
        self.covar_module_primary.base_kernel.lengthscale = init_lengthscale_primary

        # Secondary kernel
        self.covar_module_secondary = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=nu, ard_num_dims=secondary_dim)
        )
        self.covar_module_secondary.base_kernel.lengthscale = init_lengthscale_secondary

    def forward(self, x):
        mean_x = self.mean_module(x)

        # Split input
        x_primary = x[:, :self.primary_dim]
        x_secondary = x[:, self.primary_dim:self.primary_dim + self.secondary_dim]

        # Get kernels
        covar_primary = self.covar_module_primary(x_primary)
        covar_secondary = self.covar_module_secondary(x_secondary)

        covar_total = covar_primary + covar_secondary

        return gpytorch.distributions.MultivariateNormal(mean_x, covar_total)
    

class AdditiveGPModel_analyte_gpytorch(gpytorch.models.ExactGP): # this includes analyte information
    def __init__(self, train_x, train_y, likelihood,
                 primary_dim, secondary_dim,
                 nu=2.5,
                 init_lengthscale_primary=5,
                 init_lengthscale_secondary=5):
        
        super(AdditiveGPModel_gpytorch, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()

        self.primary_dim = primary_dim
        self.secondary_dim = secondary_dim

        self.covar_module_primary = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=nu, ard_num_dims=primary_dim)
        )
        self.covar_module_primary.base_kernel.lengthscale = init_lengthscale_primary

        # Secondary kernel
        self.covar_module_secondary = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=nu, ard_num_dims=secondary_dim)
        )
        self.covar_module_secondary.base_kernel.lengthscale = init_lengthscale_secondary

    def forward(self, x):
        mean_x = self.mean_module(x)

        # Split input
        x_primary = x[:, :self.primary_dim]
        x_secondary = x[:, self.primary_dim:self.primary_dim + self.secondary_dim]

        # Get kernels
        covar_primary = self.covar_module_primary(x_primary)
        covar_secondary = self.covar_module_secondary(x_secondary)

        covar_total = covar_primary + covar_secondary

        return gpytorch.distributions.MultivariateNormal(mean_x, covar_total)