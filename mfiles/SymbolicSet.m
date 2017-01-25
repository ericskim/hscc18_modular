classdef SymbolicSet < handle 
% two matlab commands to load grid points from a file stored by
% write_to_file(const SymbolicSet& set, const BDD& bdd, const std::string& filename, char mode='B')
% 
% Please have a look in ./manual/manual.pdf for usage details.
%
  properties (SetAccess=private)
    filename  % name of file which contains the SymbolicSet
    h         % C++ object handles
    dim       % dimension 
  end
  methods 
    function obj=SymbolicSet(filename)
    % the constructor opens the file 
      if(isstr(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end

      [h dim]=mexSymbolicSet('init',filename);
      obj.filename=filename;
      obj.h = h;
      obj.dim = dim;
      if(isempty(h))
        error(['SymbolicSet: could not read SymbolicSet from ', filename])
      end
    end
    function disp(obj)
      disp(['Matlab object to access the points cointained in the SymbolicSet stored in ', obj.filename])
      disp(' ')
    end
    function delete(obj)
      if(~isempty(obj.h))
        mexSymbolicSet('delete', obj.h);
      end
    end

    function points=restriction(obj,x,varargin)
      x=x(:);
      n=length(x);
      if(n>=obj.dim)
        error('SymbolicSet: vector equals or exceeds dimension');
      end
      if(isempty(varargin))
        dim=0:1:x-1;
        dim=dim(:);
      else 
        dim=varargin{1}(:);
        if(length(dim)~=n)
          error('SymbolicSet: vector and project dimensions must be of the same length');
        end
        if(max(dim)>obj.dim)
          error('SymbolicSet: project dimensions exceed dimension');
        end
        dim = dim-1;
      end
      points=mexSymbolicSet('restriction',obj.h,x,dim);
      if(isempty(points))
          error(['SymbolicSet: no grid points found for ',num2str(x)]);
      end
    end

    function points=points(obj,varargin)
    % read the grid points from file
    % optional projection indices
    % get.point('project',[1 3]) loads the grid points in the symbolic set
    % projected onto the dimensions [1 3]
      if(isempty(varargin))
        dim=0:1:obj.dim-1;
        dim=dim(:);
      else 
        dim=varargin{1}(:);
        if(length(dim)>=obj.dim)
          error('SymbolicSet: vector of project dimensions exceeds dimension');
        end
        if(max(dim)>obj.dim)
          error('SymbolicSet: project dimensions exceed dimension');
        end
        dim = dim-1;
      end
      points=mexSymbolicSet('gridpoints',obj.h,dim);
    end
  end
end
