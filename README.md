# normApp

## Description
 L'application normApp est une application Rshiny destiné à corriger les effets batch des fichiers à partir de billes de QC. L'idée est de réajuster les intensités de fluorescence de fichiers à partir de celles des billes de fichier de référence. 

## App sections
- **Upload**: Selectionner le nombre de groupes, un groupe = un groupe de fichier qui proviennent du meme batch. Importer vos fichiers : le groupe 1 correspond au groupe de référence. Vous pouvez importer : soit des fichiers contenant les cellules + les billes au sein du meme fichier soit séparés.
  
- **Preprocessing**: Apply compensation and transformation to the data.
- **Beads extraction**: 
  - ***Preprocessing***: Cette partie permet de compenser et de transformer vos données. Je vous conseille de ne pas compenser ni transformer. 
  - ***Mini-batch K-means***: Clusterisation de vos fichiers à l'aide de l'algorithme mini batch kmeans
  - **Separate beads from cells**: Vous devez choisir le ou les clusters qui correspondent à vos billes (généralement SSC-A high) 
- **Bead QC**: Nettoyer les billes extraites : supprimer les doublets, et choisir un type de bille utilisé pour la normalisation.
  
- **Normalization**: Appliquer la normalisation et visualiser les résultats

## Installation
cellAnnotationApp can be installed and run using Docker. You can install it in your own laptop.

Follow the procedure :

**1.** Install docker
[ttps://docs.docker.com/desktop/install/](https://docs.docker.com/engine/install/)

**2.** Clone the cellAnnotationApp repertory
   
  ```sh
git clone https://github.com/cipheLab/cellAnnotationApp.git
  ```

**3.** Go to the repertory (replace 'path/To/Repertory' by the path to cellAnnotationApp) : 
  ```sh
cd path/To/Repertory/cellAnnotationApp
  ```
**4.** Build docker image :
  ```sh
docker build -t cell_annotation_app
  ```
**5.** Launch the docker image :
  ```sh
docker run -p 3838:3838 cellAnnotationApp
  ```
**6.** Open an internet browser and type :  "http://127.0.0.1:3838". CellAnnotationApp will appear in the bowser.
   

   ## Usage

### Upload and View Data
 
 Click on the "Upload..." button to select and upload one or multiple FCS files.

### Preprocessing
Go to the **Preprocessing** tab.
Select the markers to transform and the type of transformation (logicle or arcsinh).
Click "Apply" to perform the preprocessing.

### Annotation

#### XGBoost
Navigate to the **XGBoost** tab.
Load or choose an XGBoost model.
Select markers used for prediction.
Click "Annotate Selected Files" to perform the annotation.

#### Scyan
Go to the **Scyan** tab.
Load or build a Scyan knowledge table.
3. Select markers for prediction.
4. Click "Run Scyan Annotation" to perform the annotation.

#### Scaffold
1. Navigate to the **Scaffold** tab.
2. Perform clustering or upload already clustered files.
3. Build or load a Scaffold map.
4. Click "Run Scaffold Annotation" to perform the annotation.

### Results Visualization and Download
1. Go to the **Results** tab to view annotation statistics and visualizations.
2. Click "Download ALL CSV results" to export the results.

Download your FCS file(s) with a new column that contain the annotation and also the dimension reduction.

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request or open an issue.


## Contact
If you have any questions, feel free to reach out to [maelle-marine.monier@inserm.fr].

