function wtresult = Superlet_Parfor(input, wavelets, F, Npoints, Nbuffers, padding, order_frac, order_int)

% the output scalogram
wtresult = zeros(numel(F), Npoints);

% convenience indexers for the zero-padded buffer
bufbegin    = padding + 1;
bufend      = padding + Npoints;

% the zero-padded buffer
buffer = zeros(Npoints + 2 * padding, Nbuffers);
for i_buf = 1 : Nbuffers
    buffer(bufbegin : bufend, i_buf) = input(i_buf, :);
end

% loop over the input buffers
for i_buf = 1 : Nbuffers
    parfor i_freq = 1 : numel(F)
        % pooling buffer, starts with 1 because we're doing geometric mean
        temp = ones(1, Npoints);

        % fill the central part of the buffer with input data
        % buffer(bufbegin : bufend) = input(i_buf, :);

        % get the number of integer wavelets
        n_wavelets = floor(order_frac(i_freq));

        % compute the convolution of the buffer with each wavelet in the
        % current set (integer wavelets)
        for i_ord = 1 : n_wavelets
            % restricted convolution (input size == output size)
            tempcx = conv(buffer(:,i_buf), wavelets{i_freq, i_ord}, 'same');

            % accumulate the magnitude (times 2 to get the full spectral
            % energy), pool with exponent = 1
            temp = temp .* (2 .* abs(tempcx(bufbegin : bufend)) .^ 2)';
        end

        % handle fractional exponent
        if (is_fractional(order_frac(i_freq)) && ...
            ~isempty(wavelets{i_freq, order_int(i_freq)}))
            % set the order index
            i_ord2 = order_int(i_freq);

            % the exponent is the fractional remainder
            exponent = order_frac(i_freq) - fix(order_frac(i_freq));

             % restricted convolution (input size == output size)
            tempcx = conv(buffer(:,i_buf), wavelets{i_freq, i_ord2}, 'same');

            % accumulate the magnitude (times 2 to get the full spectral
            % energy), pool with exponent = 1
            temp = temp .* ((2 .* abs(tempcx(bufbegin : bufend)) .^ 2)') .^ exponent;
        end

        % compute the order of the geometric mean
        root = 1 / order_frac(i_freq);
        temp = temp .^ root;

        % accumulate the current FOI to the result spectrum
        wtresult(i_freq, :) = wtresult(i_freq, :) + temp;
    end
end

% scale the output by the number of input buffers
wtresult = wtresult ./ Nbuffers;

end

% tell me if a number is an integer or a fractional
function res = is_fractional(x)
    res = fix(x) ~= x;
end