# auau200GeVQvector
  
###Code Authors:  
	Guannan Xie (guannanxie@lbl.gov)  
	Hao Qiu (hqiu@lbl.gov)  
	Mustafa Mustafa (mmustafa@lbl.gov)  

- - -

###How to run this code:  
```bash
# To run this code, make sure you already have centrality information,
StRefMultCorr* grefmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr();
StEventPlane*  eventPlaneMaker = new StEventPlane("eventPlaneMaker",picoDstMaker,grefmultCorrUtil);
```
###How to get information form this maker:
```bash
 #include "StRoot/StEventPlane/StEventPlane.h"
 StEventPlane*  eventPlaneMaker;
 #meventPlane(eventPlaneMaker);
 #cout<<meventPlane->getCentrality()<<endl;  
 #cout<<meventPlane->getEventPlane()<<endl;
 #cout<<meventPlane->getResolutionRandom()<<endl;
 ```
