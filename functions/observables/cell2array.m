function output =  cell2array(input, row)

[~, col] = size( input ) ;

[ri, ci] = size( input{row} ) ;

if ri > 1

    output = zeros(ri, ci, col);
        for j=1:col
            output(:,:,j) = input{row, j};
        end
else

    output = zeros(ci, col);
    for j=1:col
        output(:,j) = input{row, j};
    end

end

end