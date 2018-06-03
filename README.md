# image-analysis
Examples illustrating basic stages of analyzing calcium imaging data

### Requirements
We will use the `TIFFStack` class to load image stacks, which can be found on GitHub: 
https://github.com/DylanMuir/TIFFStack

### Contents

* RoiMaker - A simple class that makes a GUI to interactively place elliptical ROIs on an image
* dftregister - Calculates x/y frame translation to remove motion artifacts
* shiftframe - Fast method for applying frame translations
* image_analysis_pipeline - example basic image analysis pipeline