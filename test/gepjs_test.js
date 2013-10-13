'use strict';

var gepjs = require('../lib/gepjs.js');
var fs=require('fs');

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


exports['bettingodds'] = function(test) {
     var gep = new gepjs.GepStruct(paramOptions,options,paramValuesMap);     
      console.log(gep.processData.fittestValue); 
      console.log(gep.processData.fittest_ind); 
      
    test.ok(gep.processData.fittestValue>0);    
    test.done(); 
};


  


/*
  ======== A Handy Little Nodeunit Reference ========
  https://github.com/caolan/nodeunit

  Test methods:
    test.expect(numAssertions)
    test.done()
  Test assertions:
    test.ok(value, [message])
    test.equal(actual, expected, [message])
    test.notEqual(actual, expected, [message])
    test.deepEqual(actual, expected, [message])
    test.notDeepEqual(actual, expected, [message])
    test.strictEqual(actual, expected, [message])
    test.notStrictEqual(actual, expected, [message])
    test.throws(block, [error], [message])
    test.doesNotThrow(block, [error], [message])
    test.ifError(value)
*/

