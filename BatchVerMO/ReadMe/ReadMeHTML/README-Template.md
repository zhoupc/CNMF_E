# CNMF-E (BatchVer)

CNMF-E is powerful in extracting neurons from calcium imaging data captured by endoscope, especially in longer data. However, due to the concern of photo-bleaching and other experimental limits and designs, usually our data sets consist of large amount of small(short) data. One can run CNMF-E on all the small data individually yet if you would like to track the same nuerons over days, some extra analysis would be needed. This version helps deal with this situation. Most inspirations of this version derive from ideas implemented in the normal CNMF-E. It has good scalability and manageable memory use. Please read CNMF-E principle before using this version as the documentation assumes the terminology used in normal CNMF-E. See section "Description" for method descirption, and see "Getting Started" to use this version of cnmf-e.

## Descirption
The package uses normal CNMF-E on two different data sets with two different "depth":

 (1) full CNMF-E on sample data---"initiation mode"---to get a representative A of this data set; 
 
 (2) simplified CNMF-E on all data individually, regressing background subtracted data onto A (representative A found in step 1) and C---"massive mode". 
 
 The following figure shows the general idea of the batched version of CNMF-E.
 
![Figure1, General Description of Batch Version CNMF-E](/Users/gushijie/Documents/MATLAB/CaImaging/cnmf-e/BatchVer/ReadMe/General.png)

 More detailed decription of this version is shown in Figure2 below.
![Figure2, General Description of Batch Version CNMF-E](/Users/gushijie/Documents/MATLAB/CaImaging/cnmf-e/BatchVer/ReadMe/Detail.png)
 
### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

