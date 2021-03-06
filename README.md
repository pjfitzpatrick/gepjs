# gepjs

A Javascript Implementation of Gene Expression Programming.

## Getting Started
Install the module with: `npm install gepjs`

```javascript


## Documentation

var fs=require('fs');
//sample file of soccer odds data (in the test directory)
var fileStr=fs.readFileSync("Odds.json", 'utf8');
var paramValuesMap = JSON.parse(fileStr)

var options={'poplSize':100,
'numGenerations':100,
'diffFn':function(a,b) {return a-b;},
'mutationRate':0.3,
'crossOverGeneRate':0.2,
'crossOver1Rate':0.4,
'crossOver2Rate':0.4,
'transIsRate':0.2,
'transRisRate':0.2,
'transGeneRate':0.2,
'keepFittestN':5,
'numGenes':2,
'headSize':4,
'transpoLengthList':[1,2]};

var paramOptions={
resultKey:"FT_Home",
terminalMap:{"CS_00":"","CS_10":""},
funcKeyList:["+2","-2","*2","/2","max2","neg","inv","pow2","pow0.5"],
floatValue:[],
intValue:[]
};



var gep = new gepjs.GepStruct(paramOptions,options,paramValuesMap);     
console.log(gep.processData.fittestValue); 

console.log(gep.processData.fittest_ind); 


## License
Copyright (c) 2013 PJ Fitzpatrick  
Licensed under the MIT license.
