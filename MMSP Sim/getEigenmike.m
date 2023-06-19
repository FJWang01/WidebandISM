
function hom = getEigenmike()

% hom = getEigenmike()
% hom.gain = Calibrated gains for the Eigenmic in ANU AASP Lab
% hom.cart = mic positions in Cartesian coordinates
% hom.dist = Distance between mics
    
hom.rigid = 1;
hom.r = 0.042; % radius of the Eigenmike, in metres
hom.el = [1.2043; 1.5708; 1.9373; 1.5708; 0.5585; 0.9599; 1.5708; 2.1817; 2.5831; 2.1817; 1.5708; 0.9599; 0.3665; 1.0123; 2.1118; 2.7751; 1.2043; 1.5708; 1.9373; 1.5708; 0.5585; 0.9599; 1.5708; 2.1817; 2.5831; 2.1817; 1.5708; 0.9599; 0.3665; 1.0123; 2.1293; 2.7751];
hom.az = [0; 0.5585; 0; 5.7247; 0; 0.7854; 1.2043; 0.7854; 0; 5.4978; 5.0789; 5.4978; 1.5882; 1.5708; 1.5708; 1.5533; 3.1416; 3.7001; 3.1416; 2.5831; 3.1416; 3.927; 4.3459; 3.927; 3.1416; 2.3562; 1.9373; 2.3562; 4.6949; 4.7124; 4.7124; 4.7298];
%hom.cart = SHTools.s2c(hom.el, hom.az, hom.r);
%hom.dist = squareform(pdist(hom.cart));
%hom.w = 1;

% Calibrated gain of the Eigenmike at ANU AASP lab
%hom.gain = [1;0.841783705978144;0.824422770634782;1.07805634428050;0.819368701472337;0.934155064445381;0.942818510614020;0.876545989186322;1.10264739127246;0.848483639170549;1.04461665900460;0.936951042253091;0.922256134152427;1.04479150187472;0.930715983737848;1.27539079246131;1.13043435730242;1.16837435166048;0.862191048272407;0.861692088104175;1.00844381347731;0.754351026481360;0.973011316700164;0.781515901903590;0.787410743970355;0.735444383147764;0.703536416735698;0.738175096862642;0.747236005060125;0.758104189173352;0.742033871663043;1.02480445158795];
end