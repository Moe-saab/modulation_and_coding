 function v_nodes = hardDecoding(r,H, max_iterations)

[K,N] = size(H);

v_nodes=r; %initially variable nodes take values from the received vector
syndrome = v_nodes*H';
iterations=0;
while(iterations<max_iterations && norm(mod(syndrome,2))~=0)
    % STEP 1 : send v_nodes to c_nodes
    v_to_c=zeros(K,N);
    for i = 1:N
       index = find(H(:,i)); %find the indices of non zero entries at column i
       %check node f(j) connected to the variable node c(i) if element H(ji) = 1
       v_to_c(index,i)=v_nodes(i);
    end
    %in the end v_to_c is a matrix in which each row represent a check node
    % composed of received c(i)'s [0 or 1] at the connected nodes and 0's elsewhere

    % STEP 2 : the message from the c_node to the v_nodes is calculated 
    %          from the previously sent message
    c_to_v=zeros(K,N);
    for k=1:K
        index = find(H(k,:)); %find the indices of non zero entries at row k
        for i = 1:length(index)
            new_index=index;
            new_index(i)=[];% exclude the i'th element from index vector
            c_to_v(k,index(i))= mod(sum(v_to_c(k,new_index)),2);
        end
    end
    %in the end c_to_v is a matrix in which each row represent a check node
    % sent data to variable nodes
    
    % STEP 3 : majority voting
    for c=1:N
        index = find(H(:,c)); %find the indices of non zero entries at column c
        vote=(sum(c_to_v(index,c))+v_nodes(c))/(length(index)+1);
        if(vote>0.5)
            v_nodes(c)=1;
        else
            v_nodes(c)=0;
        end   
    end

    syndrome=v_nodes*H.';
    iterations = iterations +1;
end
end