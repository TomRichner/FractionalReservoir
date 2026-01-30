function varargout = gpuReservoirOps(operation, varargin)
%GPURESERVOIROPS GPU-accelerated operations for reservoir computing
%
% This function provides GPU support for matrix operations while keeping
% ODE integration on CPU (hybrid approach).
%
% Usage:
%   use_gpu = gpuReservoirOps('check');
%   W_gpu = gpuReservoirOps('convert', W);
%   result = gpuReservoirOps('multiply', W, states);
%   W_cpu = gpuReservoirOps('gather', W_gpu);
%
% Operations:
%   'check' - Check if GPU is available
%   'convert' - Convert matrix to gpuArray
%   'multiply' - Matrix multiplication (W * states)
%   'gather' - Transfer gpuArray back to CPU
%   'clear' - Clear GPU memory

    persistent gpu_available;
    persistent gpu_initialized;
    
    % Initialize GPU check once
    if isempty(gpu_available)
        try
            if canUseGPU
                gpu_available = true;
                gpu_initialized = true;
            else
                gpu_available = false;
                gpu_initialized = false;
            end
        catch
            gpu_available = false;
            gpu_initialized = false;
        end
    end
    
    switch lower(operation)
        case 'check'
            varargout{1} = gpu_available;
            
        case 'convert'
            % Convert matrix to gpuArray if GPU available
            if gpu_available
                try
                    varargout{1} = gpuArray(varargin{1});
                catch
                    % Fallback to CPU if conversion fails
                    varargout{1} = varargin{1};
                end
            else
                varargout{1} = varargin{1};
            end
            
        case 'multiply'
            % Matrix multiplication: result = A * B
            A = varargin{1};
            B = varargin{2};
            
            if gpu_available
                try
                    % Convert to GPU if not already
                    if ~isa(A, 'gpuArray')
                        A = gpuArray(A);
                    end
                    if ~isa(B, 'gpuArray')
                        B = gpuArray(B);
                    end
                    
                    % Perform multiplication on GPU
                    result = A * B;
                    
                    % Keep result on GPU (caller can gather if needed)
                    varargout{1} = result;
                catch
                    % Fallback to CPU if GPU operation fails
                    varargout{1} = gather(A) * gather(B);
                end
            else
                % CPU fallback
                varargout{1} = A * B;
            end
            
        case 'gather'
            % Transfer gpuArray back to CPU
            if isa(varargin{1}, 'gpuArray')
                varargout{1} = gather(varargin{1});
            else
                varargout{1} = varargin{1};
            end
            
        case 'clear'
            % Clear GPU memory
            if gpu_available
                try
                    reset(gpuDevice);
                catch
                    % Ignore errors
                end
            end
            
        otherwise
            error('gpuReservoirOps:UnknownOperation', ...
                'Unknown operation: %s', operation);
    end
end

function available = canUseGPU
    %CANUSEGPU Check if GPU computing is available
    try
        available = parallel.gpu.GPUDevice.isAvailable;
    catch
        available = false;
    end
end

