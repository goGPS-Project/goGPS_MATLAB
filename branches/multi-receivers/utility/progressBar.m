function f = progressBar(nMax,str)
% ProgressBar    Ascii progress bar.
%  progbar = progressBar(nmax,str) creates a single-line 0...10 progress bar and returns a
%  pointer to a function handle which can then be called to update it and show in the command window.
%  Both elapsed and remaining time are shown.
%{
  to update, call progbar(currentstep)

  example:
     steps = 500;
     progbar = progressBar(steps,'Rendering');
     for tmp = 1:steps
       progbar(tmp);
       pause(.01)
     end
  
  Revisions:
  version 1.0, 2012-12-08 by Alan Robinson, heavily adapted code from:

%  by david szotten 2008
%   merged with utility by us / cssm by: dn 2008
%  2008-09-16  DN Added elapsed time and estimated time left
%}


lastPercentileWritten = 0;
pstrlen = 0;
tstrlen = 0;


fprintf('%s ',str); % print out lable
hllen = length([str ' ']);

t=datenum(clock);

f = @updateBar;
    function updateBar(nCurrent)
        
        %what percentile are we up to
        currentPercentile = round(30*nCurrent/nMax);
        
        fprintf('%s',repmat(char(8),1,tstrlen)); % remove time string
        
        % compute time info
        ttn = datenum(clock)-t;
        tt  = datevec(ttn);
        dtt = ttn/nCurrent;
        ttleft = datevec(dtt*(nMax-nCurrent));
        tstr = sprintf(' time: %d:%d:%d, left: %d:%d:%d',tt(4),tt(5),round(tt(6)),ttleft(4),ttleft(5),round(ttleft(6)));
        tstrlen = length(tstr);

%         %have we passed another percentile?
%         if (currentPercentile > lastPercentileWritten )
% 
%             %we may have increased by several percentiles,
%             %so keep writing until we catch up
%             percentileToWrite = lastPercentileWritten + 1;
%             while(percentileToWrite <= currentPercentile)
% 
%                 %for every 10th, use a '+' instead
%                 if( mod(percentileToWrite,3)==0 )
%                     fprintf('%i',floor(percentileToWrite/3));
%                     
%                 else
%                     fprintf('%s','.');
%                 end
%                 percentileToWrite = percentileToWrite + 1;
%                 pstrlen = pstrlen + 1;
%             end
% 
%             %update status
%             lastPercentileWritten = currentPercentile;
% 
%         end
            
        fprintf('%s',tstr); % write time string
        
        
    end

end