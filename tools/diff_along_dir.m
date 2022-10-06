function y = diff_along_dir( x, dir, bc_id )
%DIFF_ALONG_DIR Returns directional gradient along a direction of a 3d image
% For a matrix, x, this function returns the gradient of 'x' along any of
% the 26 directions.
%
% Directions (row, col, stack):     1 = ( 0,  1,  0) % xy plane
%                                   2 = (-1,  1,  0)
%                                   3 = (-1,  0,  0) 
%                                   4 = (-1, -1,  0)
%                                   -1 = ( 0, -1,  0)
%                                   -2 = ( 1, -1,  0)
%                                   ...
%
% ---- INPUTS ----
% x       = the input image matrix
% dir     = the direction id number
% bc_id   = boundary condition identifier
%           1 = periodic (default)
%           2 = zero-pad
%           3 = replicate
%--------------------------------------------------------------------------
    if ~exist('bc_id','var')
        bc_id = 1;
    end
    switch bc_id
        case 1
            shift_op = @(x,shift)(circshift(x,shift));
        case 2
            shift_op = @(x,shift)(zeroshift(x,shift));
        otherwise
            error('unsupported boundary condition identifier');
    end
    
    switch dir
    case 1
        y = shift_op(x,[0,-1,0]) - x;
    case 2
        y = shift_op(x,[1,-1,0]) - x;
    case 3
        y = shift_op(x,[1,0,0]) - x;
    case 4
        y = shift_op(x,[1,1,0]) - x;
    case 5
        y = shift_op(x,[0,-1,-1]) - x;
    case 6
        y = shift_op(x,[1,-1,-1]) - x;
    case 7
        y = shift_op(x,[1,0,-1]) - x;
    case 8
        y = shift_op(x,[1,1,-1]) - x;
    case 9
        y = shift_op(x,[0,1,-1]) - x;
    case 10
        y = shift_op(x,[0,0,-1]) - x;
    case 11
        y = shift_op(x,[1,-1,1]) - x;
    case 12
        y = shift_op(x,[1,0,1]) - x;
    case 13
        y = shift_op(x,[1,1,1]) - x;
    case -1
        y = shift_op(x,[0,1,0]) - x;
    case -2
        y = shift_op(x,[-1,1,0]) - x;
    case -3
        y = shift_op(x,[-1,0,0]) - x;
    case -4
        y = shift_op(x,[-1,-1,0]) - x;
    case -5
        y = shift_op(x,[0,1,1]) - x;
    case -6
        y = shift_op(x,[-1,1,1]) - x;
    case -7
        y = shift_op(x,[-1,0,1]) - x;
    case -8
        y = shift_op(x,[-1,-1,1]) - x;
    case -9
        y = shift_op(x,[0,-1,1]) - x;
    case -10
        y = shift_op(x,[0,0,1]) - x;
    case -11
        y = shift_op(x,[-1,1,-1]) - x;
    case -12
        y = shift_op(x,[-1,0,-1]) - x;
    case -13
        y = shift_op(x,[-1,-1,-1]) - x;
    otherwise
        error('Invalid direction!')
    end

end

