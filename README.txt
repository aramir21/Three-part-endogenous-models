Files to perform simulations and applications of three-part endogenous models
1) SimulationV3Bnew.R --> (Baseline)
	DGP: (THREE PRODUCTS, J=3) Three part incidental truncation dependence between access, selection and outcome equations with diferent regressors in the three stages, although same regressors for each variable within each stage
	Estimation as DGP
2) SimulationV2Cnew.R --> (Exogenous access)
	DGP: (THREE PRODUCTS, J=3) Three part incidental truncation dependence between access, selection and outcome equations with diferent regressors in the three stages, although same regressors for each variable within each stage
	Estimation: Asume access is independent of selection and outcome equations
3) SimulationV3Cnew.R --> (NO exclusion restrictions)
	DGP: (THREE PRODUCTS, J=3) Three part incidental truncation dependence between access, selection and outcome equations with regressor in access, selections and outcomes equal to W
	Estimation as DGP
4) SimulationV1B.R -->
	Univariate setting
	DGP: Three part incidental truncation dependence between access, selection and outcome equations with diferent regressors in the three stages
	Estimation as DGP
5) SimulationV1Cnew.R --> (ONE product)
	DGP: (THREE PRODUCTS, J=3) Three part incidental truncation dependence between access, selection and outcome equations with diferent regressors in the three stages, although same regressors for each variable within each stage
	Estimation in an univariate setting
6) AppMarijuanaColV1.R --> Univariate marijuana application
7) AppSSPPCol.R --> Multivariate SSPP application
8) SimulationOneProduct3Stages.R -->
	DGP: ONE PRODUCT Three part incidental truncation dependence between access, selection and outcome equations with diferent regressors in the three stages, although same regressors for each variable within each stage
	Estimation as DGP
9) SimulationOneProduct3StagesExAccess.R -->
	DGP: ONE PRODUCT Three part incidental truncation dependence between access, selection and outcome equations with diferent regressors in the three stages, although same regressors for each variable within each stage
	Estimation assuming exogenous access
10) SimulationOneProduct3StagesNOexcRest.R -->
	DGP: ONE PRODUCT Three part incidental truncation dependence between access, selection and outcome equations with diferent regressors in the three stages, although same regressors for each variable within each stage
	Estimation without exclusion restrictions
 