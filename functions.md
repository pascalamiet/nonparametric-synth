**Script Author(s): PASCAL AMIET**

**Script Created on: 12.08.2023**

*Last Updated: 13.01.2024*

------------------------------------------------------------------------

**synthCtrl()**

*Description:* : This function computes the synthetic control group
given some input data and a given method.

*Parameters:*

-   `data`: The first input is the data matrix that is to be evaluated.
    It is important that the rows are the periods and the columns are
    the different groups/covariates.
-   `treated`: Here, we define the column which will be treated,
    i.e. the treatment group. This is either a string or the column
    number.
-   `t`: Here we indicate the treatment period. This has to correspond
    to the specific index.
-   `method`: The method defines which statistical model is to be used
    for the computations. One can choose between `lasso`, `ols`, `gam`.
-   `significance`: Indicates whether or not to calculate the
    significance level of the estimate. The default is
    `significance = TRUE`.

*Return Object:* The function returns a list with the estimation of the
synthetic control group, as well as some error statistics, significance,
and the absolute treatment effect.

------------------------------------------------------------------------

**getStats()**

*Description:* This function calculates all the necessary numbers to
compare methods.

*Parameters:*

-   `xtrain`: The training covariates.
-   `ytrain`: The training response variable.
-   `xtest`: The testing covariates.
-   `ytest`: The testing response variable.
-   `synth`: The predictions for the synthetic control group.

*Return Object:* The function returns a matrix containing MSE, RMSE,
MAE, and MAPE for both in-sample predictions and out-of-sample
predictions.

------------------------------------------------------------------------

**testStat()**

*Description:* This function computes the test statistic, i.e. the
significance level, of for the synthetic control group. It does so by
using a placebo test as suggested by the literature.

*Parameters:*

-   `data`: The first input is the data matrix that is to be evaluated.
    It is important that the rows are the periods and the columns are
    the different groups/covariates.
-   `treated`: Here, we define the column which will be treated,
    i.e. the treatment group. This is either a string or the column
    number.
-   `t`: Here we indicate the treatment period. This has to correspond
    to the specific index.
-   `method`: The method defines which statistical model is to be used
    for the computations. One can choose between `lasso`, `ols`, `gam`.
-   `plot`: If true, it will return a histogram of the estimated gaps to
    illustrate the significance. The default is that `plot = TRUE`.

*Return Object:* The function returns the significance level provided by
the placebo test.

------------------------------------------------------------------------

**treatEff()**

*Description:* This function calculates the treatment effect given the
treatment and the (synthetic) control group.

*Parameters:*

-   `treatment`: Time series of the treatment group.
-   `control`: Time series of the (synthetic) control group.
-   `t`: Treatment period.

*Return Object:* The function returns the treatment effect.

------------------------------------------------------------------------

**baggedTree()**

*Description:* This function fits a bagged tree model on the data.

*Parameters:*

-   `data`: The data to be fitted.
-   `y`: The variable of interested.
-   `B`: The number of trees to be fitted, by default `B = 20`.
-   `contr`: Controls for the splits when estimating the single trees,
    given as the standard `rpart.control(...)` element from the `rpart`
    package.

*Return Object:* The function returns the fitted bagged tree model.

------------------------------------------------------------------------

**predict.baggedTree()**

*Description:* This function is used for prediction of the bagged tree
model.

*Parameters:*

-   `model`: The bagged tree model that one gets as an output from the
    `baggedTree()` function.
-   `newdata`: The data (covariates) used to predict.

*Return Object:* The function returns a vector of predicted variables.

------------------------------------------------------------------------

**boostedTree()**

*Description:* This function fits a boosted tree model on the data.

*Parameters:*

-   `data`: The data to be fitted.
-   `y`: The variable of interest.
-   `B`: The number of trees to be fitted, by default `B = 20`.
-   `rate`: The learning rate, by default `rate = 0.1`.
-   `info`: Whether or not information should be displayed during
    training, by default `infor = TRUE`.
-   `contr`: Controls for the splits when estimating the single trees,
    given as the standard `rpart.control(...)` element from the `rpart`
    package.

*Return Object:* The function returns the fitted boosted tree model.

------------------------------------------------------------------------

**predict.boostedTree()**

*Description:* This function is used for prediction of the boosted tree
model.

*Parameters:*

-   `model`: The bagged tree model that one gets as an output from the
    `boostedTree()` function.
-   `newdata`: The data (covariates) used to predict.

*Return Object:* The function returns a vector of predicted variables.

------------------------------------------------------------------------

**dissimilarityMatrix()**

*Description:* This function uses distance measures to estimate the
weights of different control groups in a non-parametric way.

*Parameters:*

-   `X`: The data frame used for fitting the weights.
-   `k`: The variable of interest as a string.
-   `id`: The variable identifying the different groups as a string.
-   `id.interest`: The group we want to create a counterfactual for as a
    string.
-   `periods`: The variable name as a string of the different periods.
-   `intervention`: The period in which the intervention took place, as
    an integer.
-   `controls`: The control variables as strings in a vector.
-   `distance`: The distance measure to compute the disimilarity between
    observations. The options are
-   `method`: The method of selecting a single weight per control group.
    This can be either `'Aggregate'` or `'Norm'`, where aggregating
    means taking the mean. By default, `method = 'Norm'`.
-   `norm.wgt`: Boolean value indicating whether or not to normalize the
    estiamted weights, such that the sum of all weights is equal to 1.
    The default is set to `norm.wgt = FALSE`.
-   `true_cf`: A vector containing the true counterfactual, if that is
    known. This will allow to calculate the MSPE. This does not
    necessairily have to be provided.
-   `standardize`: A boolean value indicating whether or not to
    standardize the control variables, i.e. setting mean to 0 and
    standard deviation to 1. By default, this is not done,
    i.e. `standardize = FALSE`.

*Return Object:* The function returns a list. The list contains the
dissimilarity matrix, the estimated weights, predictions, as well as
performance measures.

------------------------------------------------------------------------

**uniKernel()**

*Description:* This function uses kernels to estimate the weights of
different control groups in a non-parametric way. It uses only
univariate kernels.

*Parameters:*

-   `X`: The data frame used for fitting the weights.
-   `k`: The variable of interest as a string.
-   `id`: The variable identifying the different groups as a string.
-   `id.interest`: The group we want to create a counterfactual for as a
    string.
-   `periods`: The variable name as a string of the different periods.
-   `intervention`: The period in which the intervention took place, as
    an integer.
-   `controls`: The control variables as strings in a vector.
-   `method`: The method of selecting a single weight per control group.
    This can be either `'Aggregate'` or `'Norm'`, where aggregating
    means taking the mean. By default, `method = 'Norm'`.
-   `h`: The bandwidth used for the kernel. If set to `NULL`, which is
    the default, then the function will choose the optimal bandwidth
    based on cross validation.
-   `kernel`: The kernel used. The default is `kernel = Gaussian`, but
    one can also choose to instead use `Epanechnikov`, `Triangular`, or
    `Tri-cube`.
-   `norm.wgt`: Boolean value indicating whether or not to normalize the
    estiamted weights, such that the sum of all weights is equal to 1.
    The default is set to `norm.wgt = FALSE`.
-   `true_cf`: A vector containing the true counterfactual, if that is
    known. This will allow to calculate the MSPE. This does not
    necessairily have to be provided.
-   `standardize`: A boolean value indicating whether or not to
    standardize the control variables, i.e. setting mean to 0 and
    standard deviation to 1. By default, this is not done,
    i.e. `standardize = FALSE`.

*Return Object:* The function returns a list. The first element of that
list is some specifications about the method. The second element are the
estimated weights.

------------------------------------------------------------------------

**multiKernel()**

*Description:* This function uses kernels to estimate the weights of
different control groups in a non-parametric way. It uses multivariate
kernels or the product of different univariate kernels.

*Parameters:*

-   `X`: The data frame used for fitting the weights.
-   `k`: The variable of interest as a string.
-   `h`: The bandwidth used for the kernel. If set to `NULL`, which is
    the default, then the function will choose the optimal bandwidth
    based on cross validation.
-   `kernel`: The kernel used. The default is `kernel = Gaussian`, but
    one can also choose to instead use `Epanechnikov`, `Triangular`, or
    `Tri-cube`.
-   `norm.wgt`: Boolean value indicating whether or not to normalize the
    estiamted weights, such that the sum of all weights is equal to 1.
    The default is set to `norm.wgt = FALSE`.

*Return Object:* The function returns a list. The first element of that
list is some specifications about the method. The second element are the
estimated weights.

------------------------------------------------------------------------

**klDiv()**

*Description:* This function allows us to easily compute the cross
entropy of two given vectors.

*Parameters:*

-   `x`: A vector of numerical values.
-   `y`: A vector of numerical values.
-   `h`: The bandwidth for bins in the discrete case.
-   `discrete`: A boolean value indicating whether or not to use the
    discrete version of the cross-entropy. By default, `discrete = T`.
-   `standardize`: A boolean value indicating whether or not to
    standardize the input vectors such that mean is zero and standard
    deviation is one. By default, `standardize = F`.

*Return Object:* The function returns a single numerical value.

------------------------------------------------------------------------

**predict.weighted()**

*Description:* This function is used to make predictions when weights
are used, such as in the kernel methods.

*Parameters:*

-   `X`: The covariates used for prediction.
-   `w`: The weights used for the prediction, i.e. the ones coming from
    the `uniKernel()` function.
-   `intercept`: Indicates the intercept value. If no value is set, the
    function assumes no intercept. The default is set to
    `intercept = NULL`.
-   `k`: The variable of interest, i.e. the one that should be
    predicted. If this is not provided, it is assumed that it is not in
    the data frame `X`.

*Return Object:* The function returns a vector of the predictions.
Additionally, if `!is.null(k)`, it will compute the mean squared error
of the predictions. In that case, the function returns a list with the
weights and the MSE.

------------------------------------------------------------------------

**optimalBandwidth()**

*Description:* This function estimates the optimal bandwidth using cross
validation and a range of possible bandwidths.

*Parameters:*

-   `x`: The covariates for which we want to estimate weights.
-   `y`: The variable of interest.
-   `kernel`: The kernel to be used, same as in `uniKernel()`.
-   `method`: The metheod to be used, same as in `uniKernel()`.
-   `intercept`: Whether or not to use an intercept, boolean variable.
-   `norm.wgt`: Whether or not to normalize the weights, boolean
    variable.

*Return Object:* The function returns a list containing the optimal
bandwidth and the different MSEs for the different bandwidths tested.

------------------------------------------------------------------------

**dgp()**

*Description:* This function is the data generating process (DGP). It
creates a time series with certain desired attributes that can be chosen
in the function calling.

*Parameters:*

-   `n`: The length of the series.
-   `form`: Specifies the form of the function as a string, can be
    either `linear`, `power`, `quadratic`, `cubic`, or `oscillating`.
    The default is `form = 'quadratic'`.
-   `alpha`: This can be used to specify the alpha of the AR process
    used to add noise. If this is not explicitly chosen, a random alpha
    will be selected.
-   `beta`: This can be used to specify the beta of the MA process used
    to add noise. If this is not explicitly chosen, the noise simply
    comes from an AR process.
-   `gamma`: This parameter squeezes the noise added to the time series.
    If not chosen, a random factor will be selected.
-   `plot`: If true, the generated time series will be plotted when the
    function is called. The default is that `plot = TRUE`.
-   `max.diff`: The maximal length of the range. The default is
    `max.diff = 10`. See the `squeeze()` function for a better overview.

*Return Object:* The function returns a list with the generated time
series, as well as a vector of the parameters used for the functional
form.

------------------------------------------------------------------------

**DGP()**

*Description:* The function generates a data set that can be used for
the simulation study to test different estimators.

*Parameters:*

-   `n`: The length of each time series.
-   `forms`: Specifies the forms of the functions as a string vector,
    can contain either `linear`, `power`, `quadratic`, `cubic`, or
    `oscillating`. The length of the vector automatically indicates the
    number of groups in the data. The first string in the vector is used
    for the treatment group.
-   `alpha`: This can be used to specify the alpha of the AR process
    used to add noise. If this is not explicitly chosen, a random alpha
    will be selected.
-   `beta`: This can be used to specify the beta of the MA process used
    to add noise. If this is not explicitly chosen, the noise simply
    comes from an AR process.
-   `gamma`: This parameter squeezes the noise added to the time series.
    If not chosen, a random factor will be selected.
-   `plot`: If true, the generated data will be plotted when the
    function is called. The default is that `plot = TRUE`.
-   `max.diff`: The maximal length of the range. The default is
    `max.diff = 10`. See the `squeeze()` function for a better overview.

*Return Object:* The function returns a list with the generated data, as
well as a vector of the parameters used for the functional form.

------------------------------------------------------------------------

**DGPcerulli()**

*Description:* The DGP from Cerulli (2020).

*Parameters:*

-   `N`: The amount of random variables generated per observation.
-   `n`: The amount of observations in the data set.
-   `n_covariates`: The number of covariates, by default,
    `n_covariates = 3`.
-   `functional.form`: The functional form when generating the data. Can
    be either the original one from Cerulli (2020), or the `New` one. By
    default, however, `functional.form = Cerulli`.
-   `intervention`: The intervention period. If not provided, the
    intervention period will be after two thirds of the observations.
-   `plot`: Whether or not to plot the data when it’s generated.

*Return Object:* The function returns a list containing the data and the
true counterfactual.

------------------------------------------------------------------------

**testModels()**

*Description:* This function tests the performance of different models
on the data that was generated before / on real-world data.

*Parameters:*

-   `x`: The data to be used to test the predictive performance of the
    different models. Generally, the output of the `DGP()` function is
    used for that.
-   `t`: This variable indicates the time of the policy intervention in
    the case of synthetic controls. Basically this states up to which
    period we use the data for training. Everything after will be used
    for calculating the test statistics of prediction accuracy.
-   `models`: A vector of strings indicating which models are to be
    fitted. The options here are `LASSO`, `Abadie`, `GAM`,
    `baggedTrees`, `rf`, `boostedTrees`, `uniKernel`, `multiKernel` and
    `LSTM`. If nothing is chosen, the default is to use all these
    models.
-   `plot`: A boolean variable indicating whether or not to plot the
    data with the fitted models. The default is that `plot = TRUE`.

*Return Object:* The function returns a matrix of performance statistics
for the different models fitted to the data.

------------------------------------------------------------------------

**simPlot()**

*Description:* This function is used to create the desired plots for
simulation studies, namely boxplots and line plots for the different
performance metrics used.

*Parameters:*

-   `x`: The data frame with the performance metrics of the different
    methods on simulated data.
-   `simulation`: The type of simulation that is performed as a string.
    This will be used for the axis labeling of the plots. E.g., when we
    simulate data with different time series lengths, we can put
    `simulation = Series Length`.
-   `metrics`: The performance metrics to be calculated as a string in a
    vector. The possible options are `MSE`, `RMSE`, and `MAE`. The
    default is to use all three of them.
-   `file`: If this is not specified, the plots will simply be showed,
    but not saved anywhere. However, one can give the plots a name here
    (i.e., a string value) and then the plots will automatically be
    saved in the Figures folder.
-   `aspect.ratio`: Gives the aspect ratio of the plot, by default,
    height is 4 and width is 12.

------------------------------------------------------------------------

**getTS()**

*Description:* This function generates a time series.

*Parameters:*

-   `alpha`: The alpha parameter of the AR-part of the model.
-   `beta`: The beta parameter of the MA-part of the model, if needed.
-   `n`: The desired length of the time series.

*Return Object:* The function returns the time series as a vector.

------------------------------------------------------------------------

**squeeze()**

*Description:* This function squeezes a time series to a desired range.

*Parameters:*

-   `ts`: The time series to be squeezed.
-   `max.diff`: The maximal length of the range. The default is
    `max.diff = 10`.

*Return Object:* The function returns the squeezed time series.

------------------------------------------------------------------------

**stack()**

*Description:* This function takes a data frame and stacks all the
columns on top of each other.

*Parameters:*

-   `x`: The data frame to stack.

*Return Object:* A data frame with one column indicating the name of the
previous column and one with the output variables.

------------------------------------------------------------------------

**unstack()**

*Description:* This function takes a data frame and unstacks it into
different columns.

*Parameters:*

-   `x`: The data frame to unstack.

*Return Object:* A data frame with multiple columns.

------------------------------------------------------------------------

**sigmoid()**

*Description:* The sigmoid function.

*Parameters:*

-   `x`: A scalar or numerical vector.

*Return Object:* A scalar or numerical vector of the same length as the
input.

------------------------------------------------------------------------

**inSilence()**

*Description:* This function is there to completely mute all possible
outputs in the console, e.g., warnings, messages, etc.

------------------------------------------------------------------------
