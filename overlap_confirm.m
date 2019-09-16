% 중복을 확인하는 함수
% Avector_norm는 확인 대상 벡터의 norm 값!, 
% TS_norm는 Avector가 중복되어있는지 확인하고픈 리스트의 norm 값!

function [temp_is_there] = overlap_confirm(Avector_norm,TS_norm)

temp = TS_norm - Avector_norm;
if min(abs(temp))>0
    temp_is_there = 0;
else
    temp_is_there = 1;
end
