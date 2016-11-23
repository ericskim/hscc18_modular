classdef UniformGrid < handle 
% some matlab commands to load points from a file stored by
% scots::UniformGrid::writeToFile() 
%
  properties (SetAccess=private)
    filename  % name of file which contains the UniformGrid
    dim       % dimension of the real space
    z         % measurement error bound
    points    % the gird points contained in the abstract set
    eta       % grid point distance parameter for the uniform grid
    first     % first grid point of the uniform grid
    nofgp     % number of grid points in each dimension 
    objhandle % handle to C++ UniformGrid instance (created in the constructor)
    project   % optional: project the grid points onto the specified dimension
  end
  methods 
    function obj=UniformGrid(filename,varargin)
    % the constructor opens the file and creates an instance of 
    % a UniformGrid with the data read from file 'filename'

      if(isstr(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end

      objhandle=mexUniformGrid('init',filename);
      obj.filename=filename;
      obj.objhandle=objhandle;
      
      if(nargin>1 && strcmp(varargin{1},'projection'))
        obj.project=varargin{2};
        obj.dim=length(varargin{2});
        %obj.points=mexUniformGrid('set',obj.objhandle,obj.project);
      elseif(nargin>1)
        error('could not read input arguments');
      else
        obj.project=0;
      end     
    end
    function delete(obj)
        mexUniformGrid('delete', obj.objhandle);
    end

    function disp(obj)
      disp(['Matlab class to access the abstract set cointained in the UniformGrid stored in ', obj.filename])
      disp(' ')
    end
    function points=get.points(obj)
      % return the points in the abstract set
      if ~isempty(obj.points)
         points=obj.points;
      else
        if(obj.project)
          points=mexUniformGrid('set',obj.objhandle,obj.project);
        else
          points=mexUniformGrid('set',obj.objhandle);
        end
        if(points==-1)
          points=[];
        end
        obj.points=points;
      end
    end
    function dim=get.dim(obj)
      % return the dimension of real space in which UniformGrid is embedded
      if ~isempty(obj.dim)
         dim=obj.dim;
      else
        dim=mexUniformGrid('dim',obj.objhandle);
        obj.dim=dim;
      end
    end
    function eta=get.eta(obj)
      % return the UniformGrid parameter eta
      if ~isempty(obj.eta)
         eta=obj.eta;
      else
        eta=mexUniformGrid('eta',obj.objhandle);
        obj.eta=eta;
      end
    end
    function z=get.z(obj)
      % return the measurement error bound z
      if ~isempty(obj.z)
         z=obj.z;
      else
        z=mexUniformGrid('zzz',obj.objhandle);
        obj.z=z;
      end
    end
    function first=get.first(obj)
      % return the first grid point
      if ~isempty(obj.first)
         first=obj.first;
      else
        first=mexUniformGrid('first',obj.objhandle);
        obj.first=first;
      end
    end
    function first=get.nofgp(obj)
      % return the number of grid points
      if ~isempty(obj.nofgp)
         nofgp=obj.nofgp;
      else
        nofgp=mexUniformGrid('nofgp',obj.objhandle);
        obj.nofgp=nofgp;
      end
    end
    function plot(obj,varargin)
    % plot grid points in the abstract set stored in the UniformGrid
    % options 
    % - 'projection', dim: grid points are projected on dimensions given in dim
    %
    % Example: plot('projection',[1 2])
    %          plots the projection of the grid points to dim one and two
    %
    % - 'plot dim',dim: specify the dimensions to be plotted
    % - 'fix at' ,val : specify the values of the remaining dimensions  
    %
    % Example: plot('plot dim',[1 3],'fix at', 3.2) 
    %          plots the dimension 1 and 3 of the grid points where the 2nd 
    %          coordinate is fixed at value 3.2
    %
      if(isempty(obj.points))
        obj.points=mexUniformGrid('set',obj.objhandle);
        if(isempty(obj.points))
        disp(['The abstract set is empty']);
        return;
        end
      end
      points=obj.points;
%      dim=mexUniformGrid('dim',obj.objhandle);; 
%      eta=mexUniformGrid('eta',obj.objhandle);; 
      dim=obj.dim;
      eta=obj.eta;
      % check arguments
      t=0;
      tidx=[];
      for i = 1:2:length(varargin)
        name = varargin{i};
        switch name
          case 'plot dim'
            if(~isnumeric(varargin{i+1}))
              error('Wrong arguments, see help UniformGrid/plot')
            end
            dim = length(varargin{i+1});
            plotidx=varargin{i+1};
            t=t+1;
            tidx=[tidx; i];
          case 'fix at'
            if(~isnumeric(varargin{i+1}))
              error('Wrong arguments, see help UniformGrid/plot')
            end
            fixval=varargin{i+1};
            t=t+1;
            tidx=[tidx; i];
          case 'projection'
            if(~isnumeric(varargin{i+1}))
              error('Wrong arguments, see help UniformGrid/plot')
            end
            dim=length(varargin{i+1});
            points=points(:,varargin{i+1});
            points=unique(points,'rows');
            tidx=[tidx; i];
          otherwise
        end
      end
      % remove 'plot dim', 'fix at' arguments from varargin 
      varargin([tidx tidx+1])=[];
      % if 'plot dim' is given, then 'fix at' must also be specified 
      % conversly, if  'fix at' is given, then 'plot dim' must also be specified 
      if(t==1)
        error('Options "plot dim" and  "fix at" must be given together')
      end
      if(t==2)
        % the number of plot dimension + entries in fix at, must add up to dim
        if((length(plotidx)+length(fixval))==dim)
          error('Options "plot dim" and  "fix at" must be given together')
        end
        % create vector of dimensions which are not plotted 
        rmdim=1:obj.dim;
        rmdim(plotidx)=[];
        % find indices which are removed 
        idx=find(all(abs(points(:,rmdim)-fixval)<=(eta(rmdim)/2),2));
        points=points(idx,plotidx);
      end
    

      switch(dim)
      case 1
        n=length(points(:,1));
        plot(points(:,1),zeros(n,1),varargin{:})
      case 2
        plot(points(:,1),points(:,2),varargin{:})
      case 3
        plot3(points(:,1),points(:,2),points(:,3),varargin{:})
      otherwise 
        disp('plot function is supported for 1, 2 and 3 dimensions only');
      end
    end
    function plotCells(obj,varargin)
    % plot cells associated with the grid points in the abstract set stored in the UniformGrid
    % options 
    % - 'projection', dim: grid points are projected on dimensions given in dim
    %
    % Example: plot('projection',[1 2])
    %          plots the projection of the grid points to dim one and two
    %
    % - 'plot dim',dim: specify the dimensions to be plotted
    % - 'fix at' ,val : specify the values of the remaining dimensions  
    %
    % Example: plot('plot dim',[1 3],'fix at', 3.2) 
    %          plots the dimension 1 and 3 of the grid points where the 2nd 
    %          coordinate is fixed at value 3.2
    %
      if(isempty(obj.points))
        obj.points=mexUniformGrid('set',obj.objhandle);
        if(isempty(obj.points))
        disp(['The abstract set is empty']);
        return;
        end
      end
      points=obj.points;
%      dim=mexUniformGrid('dim',obj.objhandle); 
%      eta=mexUniformGrid('eta',obj.objhandle); 
      dim=obj.dim;
      eta=obj.eta;
      % check arguments
      t=0;
      tidx=[];
      for i = 1:2:length(varargin)
        name = varargin{i};
        switch name
          case 'plot dim'
            if(~isnumeric(varargin{i+1}))
              error('Wrong arguments, see help UniformGrid/plot')
            end
            dim = length(varargin{i+1});
            plotidx=varargin{i+1};
            t=t+1;
            tidx=[tidx; i];
          case 'fix at'
            if(~isnumeric(varargin{i+1}))
              error('Wrong arguments, see help UniformGrid/plot')
            end
            fixval=varargin{i+1};
            t=t+1;
            tidx=[tidx; i];
          case 'projection'
            if(~isnumeric(varargin{i+1}))
              error('Wrong arguments, see help UniformGrid/plot')
            end
            dim=length(varargin{i+1});
            points=points(:,varargin{i+1});
            points=unique(points,'rows');
            tidx=[tidx; i];
          otherwise
        end
      end
      % remove 'plot dim', 'fix at' arguments from varargin 
      varargin([tidx tidx+1])=[];
      % if 'plot dim' is given, then 'fix at' must also be specified 
      % conversly, if  'fix at' is given, then 'plot dim' must also be specified 
      if(t==0)
        points=unique(points,'rows');
      end
      if(t==1)
        error('Options "plot dim" and  "fix at" must be given together')
      end
      if(t==2)
        % the number of plot dimension + entries in fix at, must add up to dim
        if((length(plotidx)+length(fixval))==dim)
          error('Options "plot dim" and  "fix at" must be given together')
        end
        % create vector of dimensions which are not plotted 
        rmdim=1:obj.dim;
        rmdim(plotidx)=[];
        % find indices which are removed 
        idx=find(all(abs(points(:,rmdim)-fixval)<=(eta(rmdim)/2),2));
        points=points(idx,plotidx);
      end
    

      switch(dim)
      case 2
        n=length(points(:,1));
        eh=obj.eta./2;
        for i=1:n
          x=points(i,1);
          xdata=x+[-1 1 1 -1]*eh(1);
          y=points(i,2);
          ydata=y+[-1 -1 1 1]*eh(2);
          v=[xdata(:) ydata(:)];
          f=[1 2 3 4];
          patch('faces',f,'vertices',v,varargin{:})
        end
      case 3
        n=length(points(:,1));
        eh=obj.eta./2;
        for i=1:n
          x=points(i,1);
          xdata=x+[-1 1 1 -1]*eh(1);
          xdata=[xdata xdata];
          y=points(i,2);
          ydata=y+[-1 -1 1 1]*eh(2);
          ydata=[ydata ydata];
          z=points(i,3);
          zdata=z+[-1 -1 -1 -1]*eh(3);
          zdata=[zdata zdata+2*eh(3)];
          v=[xdata(:) ydata(:) zdata(:)];
          f=[1 2 3 4;
             1 2 6 5;
             2 3 7 6;
             3 4 8 7;
             4 1 5 8;
             5 6 7 8];
          patch('faces',f,'vertices',v,varargin{:})
        end
      otherwise 
        disp('plot cell function is supported for 2 and 3 dimensions only');
      end
    end
  end
end
