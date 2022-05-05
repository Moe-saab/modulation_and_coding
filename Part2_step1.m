%% Group 33 : Mohamad Saab, Mehmet fatih, Muhammad tekin
    % Part1
%% First we find the generator matrix G
% First I make H in standard form parity check matrix by hand
% using row reduction method.
% H = [1     0     0     0     0     0    -1     1    -1     0
%      0     1     0     0     0     0     1    -1     0    -1
%      0     0     1     0     0     1     0     1     0     1
%      0     0     0     1     0     0     0     1     1     1
%      0     0     0     0     1    -1    -1     0    -2    -1];
H = [1     0     0     0     0     0     1     1     1     0
     0     1     0     0     0     1     0     1     0     0
     0     0     1     0     0     1     0     1     0     1
     0     0     0     1     0     0     0     1     1     1
     0     0     0     0     1     1     1     0     0     1];

G = par2genmat(H)%convert H to G
mod(G*H.',2) % to verify GHˆT=0 since if c is a codeword, 
             % then it will be orthogonal to each row of H.

%% Then we do encoding
% The encoder simply consists of modulo-2 multiplying blocks of bits with
% the generator matrix computed from the parity check matrix.
coderate = 0.5;
N = length(H); % H dimension is (N-K) x N
K = coderate * N;  
m = randi(2,1,K)-1      % m = message vector should have same length as
                        % number of rows of G which have size K x N
c = mod(m*G,2)          % c = Codeword  -> modulo-2 multiplying

%% send message m and receive r
randIndex = randi(N);
r = c;                  
r(randIndex)=~r(randIndex)   % make an error in the received signal
 
%% decode recieved signal
cdecoded = hardDecoding(r,H,10)







