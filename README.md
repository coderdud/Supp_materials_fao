# Supp_materials_fao

To access via Google Earth Engine code editor open this link:

URL: https://code.earthengine.google.com/53a4044bcdfb5c4301e421dbaccb8322

To maake pridiction using trained models:

To make predictions on a new predictor column matrix, X, use:   yfit = c.predictFcn(X) replacing 'c' with the name of the variable that is this struct, e.g. 'trainedModel'.  X must contain exactly 98 columns because this model was trained using 98 predictors. X must contain only predictor columns in exactly the same order and format as your training data. Do not include the response column or any columns you did not import into the app.  For more information, see <a href="matlab:helpview(fullfile(docroot, 'stats', 'stats.map'), 'appclassification_exportmodeltoworkspace')">How to predict using an exported model</a>.
