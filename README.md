# GPT-NEW
Image Matching Using Global Projection Transformation Correlation and the Acceleration Algorithm Using Integral Images

# Accelerated GPT by using integral images

## Introduction
- This program performes the experiment of image matching using standard GPT matching and its variants 
- The acceleration algorithm using integral images significantly shorten the calcuation time no only in matching processes but also the template construction processes

## Acknowledgement
- The author want to express his gratitude to all the co-authors of the research
- They are Prof. Yukihiko Yamashita and Prof. Toru Wakahara

## Usage

### Compile and Execute
```
make
./execGpt [matching method] [disttype]
```
switch match method
 1: GAT matching, the conventional GAT matching
 2: GPT matching, the conventional GPT matching by alternative calculation
 3: NGAT matching, the GAT matching with norm normalization
 4: NGPT matching, the GPT matching with norm normalization
 5: SGPT matching, the conventional enhanced GPT matching
 6: NSGPT matching, the enhanced GPT matching with norm normalization
 7: NSGPT-sHOG matching, the enhanced GPT matching with norm normalization
    associated with simplified HOG patterns
 9:  FNSGPT-Inte (NSGPT-sHOG with integral images)
 11: FGAT (F + matching method = acceleration algorithm)
 12: FGPT
 13: FNGAT
 14: FNGPT
 15: FSGPT
 16: FNSGPT
 17: FNSGPT-sHOG
 
switch template table type
 0: automatical selection of the type of window size
 1: for 8-quantized gradient directions
 2: for acceleration calculation of 8-quantized gradient directions
 3: for simplified HOG patterns
 4: for acceleration calculation simplified HOG patterns
 10: automatical selection of the type of window size using acceleration algorithm
 11: automatical selection of the type of window size using acceleration algorithm of integra images
 
## Reference
```
    @inproceedings{
        title={Theoretical criterion for image matching using GPT correlation},
        author={Shizhi Zhang, Toru Wakahara, Yukihiko Yamashita},
        booktitle={2016 23rd International Conference on Pattern Recognition (ICPR)},
        year={2016},
        DOI={10.1109/ICPR.2016.7899692}
    }
    @inproceedings{
        title={Image Matching Using GPT Correlation Associated with Simplified HOG Patterns},
        author={Shizhi Zhang, Toru Wakahara, Yukihiko Yamashita},
        booktitle={2017 7th International Conference on Image Processing Theory, Tools and Applications (IPTA)},
        year={2017},
    }
``` 

## License
Apache License, Version 2.0
