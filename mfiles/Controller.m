classdef Controller < handle 
% some matlab commands to acces the reachability controller stored by
% scots::ReachabilityGame::writeToFile() 
%
  properties (SetAccess=private)
    filename  % name of file which contains the Controller
    domain    % the gird points in the domain of the controller
    value     % vector of the values of each grid point in the domain
    conhandle % handle to C++ ReachabilityGame or SafetyGame instance (created in the constructor)
    sshandle  % handle to C++ UniformGrid instance (for the state space) (created in the constructor)
    ishandle  % handle to C++ UniformGrid instance (for the input space) (created in the constructor)
    project   % optional: project the domain (notice: not the input) onto the specified dimension
  end
  methods 
    function obj=Controller(filename,varargin)
    % the constructor opens the file and creates an instance of 
    % a ReachabilityGame or SafetyGame and tow UniformGrid classes 
    % with the data read from file 'filename'

      if(isstr(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end

      [h1 h2 h3]=mexController('init',filename);
      obj.filename=filename;
      obj.conhandle=h1;
      obj.sshandle=h2;
      obj.ishandle=h3;
      
      if(nargin>1 && strcmp(varargin{1},'projection'))
        obj.project=varargin{2};
        obj.domain=mexController('domain',obj.conhandle, obj.sshandle, obj.ishandle,obj.project);
      elseif(nargin>1)
        error('could not read input arguments');
      else
        obj.project=0;
      end     
    end
    function delete(obj)
        mexController('delete', obj.conhandle, obj.sshandle, obj.ishandle);
    end

    function disp(obj)
      disp(['Matlab class to access the controller computed with Game stored in ', obj.filename])
      disp(' ')
    end
    function domain=get.domain(obj)
      % return the grid points that are in the domain of the controller
      if ~isempty(obj.domain)
        domain=obj.domain;
      else
        if(obj.project)
          domain=mexController('domain',obj.conhandle, obj.sshandle, obj.ishandle,obj.project);
          obj.domain=domain;
        else
          domain=mexController('domain',obj.conhandle, obj.sshandle, obj.ishandle);
          obj.domain=domain;
        end
      end
    end
    function value=get.value(obj)
      % return the value of the value functions at the grid points in the domain of the controller
      if ~isempty(obj.value)
        value=obj.value;
      else
        value=mexController('value',obj.conhandle, obj.sshandle, obj.ishandle);
        if(value==-1)
          value=[];
        end
        obj.value=value;
      end
    end
    function u=input(obj,x,varargin)
      % return control input associated with grid point x  
      if nargin == 0 || isempty(varargin)
        u=mexController('input',obj.conhandle, obj.sshandle, obj.ishandle,x);
      elseif nargin>1 && strcmp(varargin{1},'projection')
        project=varargin{2};
        u=mexController('input',obj.conhandle, obj.sshandle, obj.ishandle,x,project);
      elseif nargin>1
        error('could not read input arguments');
      end 
      
      if(isempty(u))
        error(['State not in the domain of the controller stored in ', obj.filename])
      end
    end
  end
end
