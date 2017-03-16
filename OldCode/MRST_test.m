gravity reset on
G = cartGrid([1 1 30], [1 1 30]*meter^3);
G = computeGeometry(G);
rock = makeRock(G,0.1*darcy(),0.2); 
% porosity ~ 0.2
% permeability ~ 100millidarcy (mD or md)
% domain: [0,1]x[0,1]x[0,30] using regular 1x1x30 cartesian grid

T = computeTrans(G, rock);
bc = pside([], G, 'TOP', 100.*barsa());
%p = 100bar at the top of the column and no-flow conditions (v*n = 0)

%%%%%%%%%%%%%%%%%
mrstModule add incomp;
fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1024*kilogram/meter^3);

sol = incompTPFA(initResSol(G,0.0),G,T,fluid,'bc',bc);

plotFaces(G,1:G.faces.num, convertTo(sol.facePressure,barsa()));
set(gca, 'ZDir','reverse'),title('Pressure[bar]')
view(3), colorbar,set(gca,'DataAspect',[1 1 10])

%%%%%%%%%%%%%%%
x = initVariablesADI(zeros(3,1));
eq1 = [ 3, 2, -4]*x + 5;
eq2 = [ 1, -4, 2]*x + 1;
eq3 = [-2, -2, 4]*x - 6;
eq = cat(eq1,eq2,eq3);
u = -eq.jacf1gneq.val

%eq = div((T/mu)*.(grad(p) - g*rho*grad(z)));

%%%%%%%%%%%%%%%%%%%
%Homogenious model;
G = cartGrid([10 10]);
rock = makeRock(G, 200*milli*darcy, 0.2);

plotCellData(G, rock.poro,'EdgeColor','none');
colorbar('horiz'); axis equal tight;

plotCellData(G, convertTo(rock.perm,milli*darcy));
colorbar('horiz'); axis equal tight; view(3);


%%%%%%%%%%%%%%%%%%

