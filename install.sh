conda create -n DeepPhos pip python=3.7.6
conda activate DeepPhos
pip install keras==2.0.0
pip install numpy
pip install pandas
pip install matplotlib
pip install scikit-learn
pip uninstall tensorflow tensorflow-estimator
pip install tensorflow-gpu==1.15
