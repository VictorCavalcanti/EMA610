close all; clear all; clc;
[name,M15]=readop4('MAA15.OP4');
[name,K15]=readop4('KAA15.OP4');
[name,PHI15] = readop4('PHIA15.OP4');
[dof15]=xreadasetA('beam15.f06');
[wn15]=xreadcycles('beam15.f06');

[name,M30]=readop4('MAA30.OP4');
[name,K30]=readop4('KAA30.OP4');
[name,PHI30] = readop4('PHIA30.OP4');
[dof30]=xreadasetA('beam30.f06');
[wn30]=xreadcycles('beam30.f06');