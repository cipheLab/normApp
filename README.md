# normApp

## Description
 The normApp application is a Rshiny application designed to correct the batch effects of files from QC beads. The idea is to readjust the fluorescence intensities of files from those of reference file beads. 

## App sections
- Upload**: Select the number of groups, one group = one group of files from the same batch. Import your files: Group 1 corresponds to the reference group. You can import: either files containing cells + beads within the same file, or separate files.
  
- Preprocessing**: Apply compensation and transformation to the data.
- Beads extraction**: 
  - ***Preprocessing***: This part allows you to compensate and transform your data. I advise you not to compensate or transform. 
  - ***Mini-batch K-means***: Cluster your files using the mini-batch kmeans algorithm.
  - **Separate beads from cells**: Select the cluster(s) that correspond to your beads (usually SSC-A high). 
- **Bead QC**: Clean the extracted beads: remove doublets, and choose a bead type used for normalization.
  
- Normalization**: Apply normalization and visualize results


## Installation

**1.** Install docker
[ttps://docs.docker.com/desktop/install/](https://docs.docker.com/engine/install/)

**2.** Clone the cellAnnotationApp repertory
   
  ```sh
git clone https://github.com/cipheLab/normApp.git
  ```

**3.** Go to the repertory (replace 'path/To/Repertory' by the path to cellAnnotationApp) : 
  ```sh
cd path/To/Repertory/normApp
  ```
**4.** Build docker image :
  ```sh
docker build -t normApp
  ```
**5.** Launch the docker image :
  ```sh
docker run -p 3838:3838 normApp
  ```
**6.** Open an internet browser and type :  "http://127.0.0.1:3838". normApp will appear in the bowser.
   

   ## Usage


## Contributing
Contributions are welcome! Please feel free to submit a Pull Request or open an issue.


## Contact
If you have any questions, feel free to reach out to [maelle-marine.monier@inserm.fr].

