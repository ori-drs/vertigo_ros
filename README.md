# vertigo_ros

How to run vertigo_ros:

For generating data with random outliers:

`roscd /vertigo_ros/vertigo/datasets`
`rosrun vertigo_ros generateDataset.py -i sphere2500/originalDataset/sphere2500.g2o -s -n 100`

where, `-s` make the loop closures switchable and `-n` defines the number of random outliers. Since newcollege dataset has its own outliers listed at the bottom of the data, it doesn't need to add random outliers.

For NewereCollege dataset, since it already consists of outliere you just need to make the constraints switchable by:

`rosrun vertigo_ros generateDataset.py -i newcollege/newcollege_long.g2o -s -n 0`

Then run:
`rosrun vertigo_ros robustISAM2-3d -i new.g2o --adaptive --relinSkip 4 --relinThresh 0.1 -o results.csv`

where, `--adaptive` selects the approach (you can choose `--linear`, `--sigmoid`). `--relinSkip` chose the relinearisation step in iSAM2. Normally, for adaptive the best result is provided by having a relinSkip between 1 to 5. --relinThresh defines the relinearisation error threshold in iSAM2. Since shape parameter (\alpha) is a global parameter, I select relinThresh between 0.1 to 0.5 for the best result and avoiding indeterminate segfault. We might need to define a separate relinThresh specifically for \alpha in the future.

Once the program finished, a file called `resut.csv` is stored in `~/catkin_ws/src/vertigo_ros/vertigo/datasets`. Then, to see the estimated trajectory, you need to open the file in matlab and plot the trajectory:

```
temp = readtable('results.csv');
ad=table2array(temp(:,2:end));
plot3(ad(:,1),ad(:,2),ad(:,3),'b')
axis equal
```

For 2d:

`rosrun vertigo_ros generateDataset.py -i manhattan/originalDataset/Olson/manhattanOlson3500.g2o -s -n 100`

`rosrun vertigo_ros robustISAM2-2d -i new.g2o --adaptive --relinSkip 5 --relinThresh 0.5`
