/*
 * gepjs
 * https://github.com/sportssmarty/gepjs
 *
 * Copyright (c) 2013 PJ Fitzpatrick
 * Licensed under the MIT license.
 */

(function () {
var gep={};

var root = this;
if (typeof module !== 'undefined' && module.exports) {
  module.exports = gep;
}
else {
  root.gep = gep;
}

gep.funcMap={"+2":{"a":2,"fn":function(inpMap) {return inpMap[0]+inpMap[1]; }},
    "-2":{"a":2,"fn":function(inpMap) {return inpMap[0]-inpMap[1]; }},
    "*2":{"a":2,"fn":function(inpMap) {return inpMap[0]*inpMap[1]; }},
    "/2":{"a":2,"fn":function(inpMap) {return inpMap[0]/inpMap[1];}},
    "max2":{"a":2,"fn":function(inpMap) {return inpMap[0]>inpMap[1]?inpMap[0]:inpMap[1];}}, 
    "neg":{"a":1,"fn":function(inpMap) {return inpMap[0]*-1;}},
    "inv":{"a":1,"fn":function(inpMap) {return 1.0/inpMap[0];}},
    "pow2":{"a":1,"fn":function(inpMap) {return Math.pow(inpMap[0],2);}},
    "pow0.5":{"a":1,"fn":function(inpMap) {return Math.pow(inpMap[0],0.5);}}, 
    //"+":{"fn":function(inpMap) {return inpMap.reduce(function(a,b){return a+b;})}}	    
    "+":{"fn":function(inpMap) {return inpMap[0]+inpMap[1]}}      
}




gep.GepStruct=function(paramOptions,options,paramValuesMap) {
  this.resultKey=paramOptions['resultKey'];    
  this.terminalMap=paramOptions['terminalMap'];    

  this.funcKeyList=paramOptions['funcKeyList']||[];
  
  this.floatValue = paramOptions['floatValue']||[]
  this.intValue = paramOptions['intValue']||[]  
  this.linkingFunc=paramOptions['linkingFunc']||"+";
  this.paramValuesMap = paramValuesMap;
  this.pop = [];
  
  this.popSize=options['poplSize']||100;
  this.numGenerations=options['numGenerations']||100;
  this.diffFn = options['diffFn']||function(a,b) {return a-b;};
  this.mutationRate=options['mutationRate']||0.3;
  this.crossGeneRate=options['crossOverGeneRate']||0.2;
  this.crossOver1Rate=options['crossOver1Rate']||0.4;
  this.crossOver2Rate=options['crossOver2Rate']||0.4;
  this.transIsRate=options['transIsRate']||0.2;
  this.transRisRate=options['transRisRate']||0.2;
  this.transGeneRate=options['transGeneRate']||0.2;
  this.keepFittestN=options['keepFittestN']||5;
  this.numGenes =options['numGenes']||2;
  this.headSize=options['headSize']||4;
  this.transpoLengthList = options['transpoLengthList']||[1,2];
  
  this.tailSize = null;     
  this.geneLength=null;  
  this.chromoLength=null;   
  this.parsAllowed=null;      
  this.fitnessCache = {};
  this.popTotalFitness=[];
  this.processStatus="started";
  this.processData = {"generationNum":0,"fittest_ind":"","fittest_value":""};  

  this.initPop();
  this.updatePop();
}
//have different fitness calculations
//only use x% of the population for fitness calculations
gep.GepStruct.prototype._calcFitness = function(inpChromo) {
    var evalChromoMap,sseCalc,paramKey,calcDiff,counter;
    counter=0;
    sseCalc=0;
    evalChromoMap = this.evalChromo(inpChromo,this.paramValuesMap);
    if (evalChromoMap=="NaN") {return "NaN";}
    for (paramKey in evalChromoMap) {      
      counter++;
      calcDiff = this.diffFn(evalChromoMap[paramKey],this.paramValuesMap[paramKey][this.resultKey]);      
      if (isNaN(calcDiff) || calcDiff==Infinity || calcDiff==-Infinity) {
	return "NaN";
      }
      sseCalc = sseCalc+Math.abs(calcDiff);      
    }
    return (counter/sseCalc);
}

gep.GepStruct.prototype._copyChromo=function(inpChromo) {
  var retChromo=[];
  var y,item;
  for (y=0,item;item=inpChromo[y];y++) {
    retChromo.push(item);
  }
  return retChromo;
}

gep.GepStruct.prototype.getRanIdx=function(inpNum) { 
  return Math.round((inpNum-1)*Math.random()) 
}

gep.GepStruct.prototype.getRanEntry=function(inpArr) {   
  return inpArr[Math.round((inpArr.length-1)*Math.random())];
}

gep.GepStruct.prototype.selectParam=function(chromoPosition) {   
  var parsAllowedIdx = this.chromoLength % this.geneLength;
  return this.getRanEntry(this.parsAllowed[parsAllowedIdx]);
}

  
gep.GepStruct.prototype.compareInds = function() {
  that = this;
  return function(a, b) { 
    return (that.getFitness(b) - that.getFitness(a));
  }
}

gep.GepStruct.prototype._getRandContantValue = function() {
  if ((this.intValue.length>0 && this.floatValue.length>0)) {
    if (Math.random() < .5) {      
      return this.intValue[0] + parseInt(Math.round(Math.random()*(this.intValue[1]-this.intValue[0])));
    }
    else {
      return this.floatValue[0] +Math.random()*(this.floatValue[1]-this.floatValue[0]);
    }    
  }
  else if (this.intValue.length>0) {
    return this.intValue[0] + parseInt(Math.round(Math.random()*(this.intValue[1]-this.intValue[0])));
  }
  else if (this.floatValue.length>0) {
    return this.floatValue[0] +Math.random()*(this.floatValue[1]-this.floatValue[0]);
  }  
}
gep.GepStruct.prototype.selectInd = function() {
var rndNum = this.popTotalFitness[this.popSize] * Math.random();
  for (var i=0;i<this.popSize;i++) {
    if ((rndNum>=this.popTotalFitness[i]) && (rndNum<=this.popTotalFitness[i+1])) {
      return this.pop[i];
    }
  }
}

gep.GepStruct.prototype.initPop = function() {
  var terminalList,constantAllList;
  var maxArity=-1;        
  for (var i=0; i< this.funcKeyList.length;i++) {    
  if (gep.funcMap[this.funcKeyList[i]]["a"] && gep.funcMap[this.funcKeyList[i]]["a"]>maxArity) {
    maxArity=gep.funcMap[this.funcKeyList[i]]["a"];
    }
  }  
  this.tailSize = this.headSize*(maxArity-1)+1;     
  this.geneLength=this.headSize+this.tailSize;  
  this.chromoLength=this.geneLength*this.numGenes;   
  this.parsAllowed=[];  
  terminalList = [];
  constantAllList=[];
  if (this.intValue.length>0 || this.floatValue.length>0) { constantAllList.push("?");} 
  
  for (var terminalKey in this.terminalMap) {terminalList.push(terminalKey);}
  for (var i = 0; i < this.geneLength; i++) {    
    if (i<this.headSize) {
      this.parsAllowed[i]=this.funcKeyList.concat(constantAllList,terminalList);
    }
    else {
      this.parsAllowed[i]=terminalList.concat(constantAllList);
    }
  }  
  i=0;
  do {    
    var chromoArr=[];
    for (var j=0;j<this.numGenes;j++) {
      for (var k=0;k<this.geneLength;k++) {
	tmpValue=this.getRanEntry(this.parsAllowed[k]);
	if (tmpValue == "?") {tmpValue=this._getRandContantValue();}	
	chromoArr[chromoArr.length]=tmpValue; 
      }
    }
    tmpFitness=this.getFitness(chromoArr);  
    if (tmpFitness!="NaN") {
      this.pop.push(chromoArr)      
    };

  } while (this.pop.length<this.popSize)
  
  thisSortFunc =this.compareInds(); 
  this.pop.sort(thisSortFunc);
  return;
}

gep.GepStruct.prototype.updatePop = function() {
  var ind1,ind2,thisSortFunc,retArray,popCounter;

  this.popTotalFitness=[0];     
  for (var i=1;i<=this.pop.length;i++) {    
    this.popTotalFitness[i]=this.popTotalFitness[i-1]+this.getFitness(this.pop[i-1]);  
  }
  tmpPop =[]
  for (var x=0; x< this.keepFittestN;x++) {
    tmpPop.push(this._copyChromo(this.pop[x]));  
  }  
  do {
  ind1=this.selectInd();
  ind2=this.selectInd();   
  
  
  ind1 = this._copyChromo(ind1);
  ind2 = this._copyChromo(ind2);    
  retArray=[];  
  if (Math.random() < this.crossOver1Rate && tmpPop.length<this.popSize-1) {
    retArray=this.crossOver1(ind1,ind2);
  }
  else if (Math.random() < this.crossOver2Rate && tmpPop.length<this.popSize-1) {
    retArray=this.crossOver1(ind1,ind2);
  }
  else if (this.numGenes>1 && Math.random() < this.crossGeneRate && tmpPop.length<this.popSize-1) {
    retArray=this.crossGene(ind1,ind2);
  } 
  else if (Math.random() < this.transIsRate) {
    retArray=[this.transIS(ind1,ind2)];
  }
  else if (Math.random() < this.transRisRate) {
    retArray=[this.transRIS(ind1,ind2)];
  }      
  else if (Math.random() < this.mutationRate) {
    retArray=[this.mutate(ind1)];
  }
  else if (Math.random() < this.transGeneRate) {
    retArray=[this.transGene(ind1)];
  }
  else {
    ind1=this._copyChromo(ind1);
    retArray=[ind1];    
  }  
  tmpFitness=this.getFitness(retArray[0]);  
  if (tmpFitness!="NaN") {tmpPop.push(retArray[0])};
  if (retArray.length==2) {
    tmpFitness=this.getFitness(retArray[1]);  
    if (tmpFitness!="NaN") {tmpPop.push(retArray[1])};
      
  }  
  } while (tmpPop.length<this.popSize)
  
  for (popCounter=0;popCounter<this.popSize;popCounter++) {
    this.pop[popCounter]=this._copyChromo(tmpPop[popCounter]);
  }
  thisSortFunc =this.compareInds(); 
  this.pop.sort(thisSortFunc);
  this.processData.generationNum=this.processData.generationNum+1;    
  this.processData.fittestValue=this.getFitness(this.pop[0]); 
  this.processData.fittest_ind=this._copyChromo(this.pop[0]); 
}
     
  //Mutate a chromo    
gep.GepStruct.prototype.mutate = function(inpChromo) {
  var siteMutate = this.getRanIdx(this.chromoLength)%this.geneLength;    
  var lookupParsIdx = siteMutate%this.geneLength;    
  tmpValue = this.getRanEntry(this.parsAllowed[lookupParsIdx]);
  if (tmpValue == "?") {tmpValue=this._getRandContantValue();}	
  inpChromo[siteMutate]=tmpValue;
  return inpChromo; 
}      
    
  //Perform a one point crossover
gep.GepStruct.prototype.crossOver1 = function(inpChromo1,inpChromo2) {
  var siteCross = this.getRanIdx(this.chromoLength);
  var retChromo1=[],retChromo2=[];
  for (var x=0;x<this.chromoLength;x++) {      
    if (x<=siteCross) {
      retChromo1[x]=inpChromo1[x]; retChromo2[x]=inpChromo2[x];
    }
    else { 
      retChromo1[x]=inpChromo2[x]; retChromo2[x]=inpChromo1[x];
    }
  }
  return [retChromo1,retChromo2]; 
}     
    
  //Perform a two point crossover
gep.GepStruct.prototype.crossOver2 = function(inpChromo1,inpChromo2) {
  var siteCross1 = this.getRanIdx(this.chromoLength);
  var siteCross2 = this.getRanIdx(this.chromoLength);      
  var retChromo1=[],retChromo2=[];
  for (var x=0;x<this.chChromoLength;x++) {     
    if (x<=siteCross1 && x<=siteCross2) {
      retChromo1[x]=inpChromo1[x];retChromo2[x]=inpChromo2[x];
    }
    else if (x<=siteCross1 || x<=siteCross2) {
      retChromo1[x]=inpChromo2[x];retChromo2[x]=inpChromo1[x];
    }
    else {retChromo1[x]=inpChromo1[x]; retChromo2[x]=inpChromo2;}
    }         
  return [retChromo1,retChromo2]; 
}
    
  //Perform a Gene Crossover
gep.GepStruct.prototype.crossGene = function(inpChromo1,inpChromo2) {
  var siteCross = this.getRanIdx(this.numGenes);  
  var retChromo1=[],retChromo2=[];
  for (var x=0;x<this.chromoLength;x++) {     
    if (x<=(siteCross*this.geneLength)) {
      retChromo1[x]=inpChromo1[x]; retChromo2[x]=inpChromo2[x];
    }
    else { 
      retChromo1[x]=inpChromo2[x]; retChromo2[x]=inpChromo1[x];
    }
  }  
  return [retChromo1,retChromo2];
}   

  
  //Perform a transRIS
gep.GepStruct.prototype.transRIS = function(transpoChromo,destChromo) {
  var thisGuess;
  var dest_loc = this.getRanIdx(this.numGenes)*this.geneLength;
  var transpoLength=this.getRanEntry(this.transpoLengthList);
  var transpo_loc_head=this.getRanIdx(this.headSize);
  var transpo_loc_guess = this.getRanIdx(this.numGenes)*this.geneLength+transpo_loc_head;  
  var transpo_loc_final=-1;  
  searchString="";
  for (var x=0;x<transpo_loc_head;x++) {    
    thisGuess=transpo_loc_guess-x;
    searchString=searchString+"/"+transpoChromo[thisGuess];
    if (gep.funcMap[transpoChromo[thisGuess]]!=undefined && transpo_loc_final==-1) {
      transpo_loc_final=thisGuess;
    }      
  }
  if (transpo_loc_final!=-1) {
    for (var x=0;x<transpoLength;x++) { 
      destChromo[dest_loc+x]=transpoChromo[transpo_loc_final+x]; 
    }
  } 
  return destChromo; 
}        
    
  //Perform a transIS
gep.GepStruct.prototype.transIS = function(transpoChromo,destChromo) {
  var dest_loc = this.getRanIdx(this.numGenes)*this.geneLength+this.getRanIdx(this.headSize);
  var transpoLength=this.getRanEntry(this.transpoLengthList);    
  //assuming that this comes from the head
  var transpo_loc_final = this.getRanIdx(this.numGenes)*this.geneLength+this.getRanIdx(this.headSize);      
  for (var x=0;x<transpoLength;x++) { 
    destChromo[dest_loc+x]=transpoChromo[transpo_loc_final+x]; 
  }
  return destChromo;
}      
    
  //Perform transGene
gep.GepStruct.prototype.transGene = function(destChromo) {
  var gene_loc = this.getRanIdx(this.numGenes-1)+1;
  var copyInp=[];
  for (var x=0; x< this.chromoLength;x++) {
  copyInp[x] = destChromo[x];
  }
  var pos_count=0;
  for (var x=0; x< this.numGenes;x++) {
    for (var y=0;y<this.geneLength;y++) {
      if (x==0) {
	destChromo[pos_count]=copyInp[pos_count+(gene_loc*this.geneLength)];
      }      
      else if (x<=gene_loc) {
	destChromo[pos_count]=copyInp[pos_count-this.geneLength];
      }
      else if (x>gene_loc) {      
      }
    pos_count=pos_count+1;
    }  
  }
  return destChromo;  
}

//calculate the parameter positions

gep.GepStruct.prototype._getGene = function(geneNum,inpChromo) {
  return inpChromo.slice(geneNum*this.geneLength,geneNum*this.geneLength+this.geneLength);
}

gep.GepStruct.prototype.getFitness=function(inpChromo) {  
  var chromoKey,thisCalcFitness;
  chromoKey = inpChromo.join("~");    
  if (this.fitnessCache[chromoKey]!=undefined) {
  return this.fitnessCache[chromoKey];
  }
  else {
    thisCalcFitness = this._calcFitness(inpChromo);
    if (thisCalcFitness=="NaN") {return "NaN";}
    this.fitnessCache[chromoKey] = thisCalcFitness;
    return thisCalcFitness;
  }  
}



gep.GepStruct.prototype.evalChromo = function(inpChromo) {
  var geneNum,geneItem,positionList,geneCalcArray,paramCalcMap,paramKey;
  geneCalcArray=[];
  paramCalcMap = {};
  for (paramKey in this.paramValuesMap) {
    geneCalcArray=[];
    for (geneNum=0;geneNum<this.numGenes;geneNum++) {
      geneItem = this._getGene(geneNum,inpChromo);
      positionList = this._getGenePositionList(geneItem);   
      thisGeneCalc=this._evalGene(geneItem,positionList,this.paramValuesMap[paramKey]);
      if (isNaN(thisGeneCalc)|| thisGeneCalc==Infinity || thisGeneCalc==-Infinity) {return "NaN";}
      geneCalcArray.push(thisGeneCalc);    
    }
  if (this.numGenes==1) {
    chromoValue = geneCalcArray[0];
  }
  else {
    chromoValue = gep.funcMap[this.linkingFunc]["fn"](geneCalcArray);
  }
  paramCalcMap[paramKey]=chromoValue;
  }
  return paramCalcMap;
}

gep.GepStruct.prototype._getGenePositionList=function(inpGene) {
var geneItem, geneItemIdx,constantCount;
var positionList = [];
var thisArity,cumArity;
  cumArity=0;
  constantCount=0;
  for (geneItemIdx=0,geneItem;geneItem=inpGene[geneItemIdx];geneItemIdx++) {
    if (this.terminalMap[geneItem]!=undefined) {
      positionList.unshift("param");
      thisArity=-1;
    }      
    else if (gep.funcMap[geneItem]!=undefined) {
      positionList.unshift("func");
      thisArity=gep.funcMap[geneItem]["a"]-1;
    }
    else {
      positionList.unshift("other");
      thisArity=-1;
    }
    cumArity=cumArity+thisArity;    
    if (cumArity<=-1) {	
	break;
    }
  }
  return positionList;
}

gep.GepStruct.prototype._evalGene=function(inpGene,positionList,paramMap) {
var positionItem, positionItemIdx,genePosition,geneItem;
var positionListLength=positionList.length;
var currentPosition,paramCounter,variableList,variableQueue;
  variableQueue=[];
  for (positionItemIdx=0,positionItem;positionItem=positionList[positionItemIdx];positionItemIdx++) {
    genePosition = positionListLength-positionItemIdx-1;
    geneItem = inpGene[genePosition];
      if (positionItem =="func") {     
	variableList=[];
	for (paramCounter=1;paramCounter<=gep.funcMap[geneItem]["a"];paramCounter++) {
	  thisVariable=variableQueue.shift();
	  variableList.push(thisVariable);
	}
	variableQueue.push(gep.funcMap[geneItem].fn(variableList));        
      }
      else if (positionItem=="param") {
	variableQueue.push(paramMap[geneItem]);
      }
      else {
	variableQueue.push(geneItem);      
      }
  }
  return variableQueue[0];
}


})();