% �ߺ��� Ȯ���ϴ� �Լ�
% Avector_norm�� Ȯ�� ��� ������ norm ��!, 
% TS_norm�� Avector�� �ߺ��Ǿ��ִ��� Ȯ���ϰ��� ����Ʈ�� norm ��!

function [temp_is_there] = overlap_confirm(Avector_norm,TS_norm)

temp = TS_norm - Avector_norm;
if min(abs(temp))>0
    temp_is_there = 0;
else
    temp_is_there = 1;
end
