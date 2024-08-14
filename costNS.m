% This function computes the real cost of extraction given the stock q of resources
function cos = costNS(Q)
global pf
cos =1000000*(0.00466226*Q^2-29.62811308*Q+80971.56768)/pf; % see explanation in mainprogramexhaust. 
end
%pf represents the price factor or price level of the resource. 
%It is used to normalize or adjust the cost calculation, ensuring that the costs are expressed in real terms rather than nominal terms. 
%This adjustment accounts for inflation or other changes in the price level over time.
%Dividing by pf adjusts the cost to reflect real terms. This normalization ensures that the cost calculation accounts for changes in the price level, 
%providing a consistent basis for comparison across different time periods or scenarios. 


