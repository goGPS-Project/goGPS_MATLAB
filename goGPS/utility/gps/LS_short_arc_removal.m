function [A, y0, b, Q, sat_track, amb_num, amb_prn_track] = LS_short_arc_removal(A, y0, b, Q, sat_track, amb_num, amb_prn_track, min_arc)

% %remove observations without pivot OR slave satellite
% rem_obs = find(sum(A(:,4:end),2));
% rem_obs = setdiff(rem_obs, ok_obs);
% if (~isempty(rem_obs))
%     A(rem_obs,:) = [];
%     y0(rem_obs) = [];
%     b(rem_obs) = [];
%     Q(rem_obs,:) = []; Q(:,rem_obs) = [];
%     sat_track(rem_obs,:) = [];
%     ok_obs = find(sum(A(:,4:end),2));
% end

%remove ambiguity unkowns with arcs shorter than given threshold
rem_amb = find(sum(A~=0,1) < min_arc);
if (~isempty(rem_amb))
    for r = 1 : length(rem_amb)
        rem_obs = find(A(:,rem_amb(r))~=0);
        A(rem_obs,:) = [];
        y0(rem_obs) = [];
        b(rem_obs) = [];
        Q(rem_obs,:) = []; Q(:,rem_obs) = [];
        sat_track(rem_obs,:) = [];
    end
    A(:,rem_amb) = [];
    amb_num = amb_num - length(rem_amb);
    amb_prn_track(rem_amb-3) = [];
    
%     %check again and remove observations without pivot OR slave satellite
%     rem_obs = find(sum(A(:,4:end),2));
%     A(rem_obs,:) = [];
%     y0(rem_obs) = [];
%     b(rem_obs) = [];
%     Q(rem_obs,:) = []; Q(:,rem_obs) = [];
%     sat_track(rem_obs,:) = [];
%     
%     rem_amb = find(sum(A~=0,1) < min_arc);
end
