function [] = im3dscroll(IM,r,tt)
%This function takes in a 3D image IM and 3-column position vector r. r can
%be empty. It allows the user to scroll through different planes, with
%points in the current plane highlighted.
%Examples: 
%im3dscroll(rand(200,200,10),[(1:200)',(1:200)',(0.55:0.05:10.50)']);caxis([0 5]);
%im3dscroll(rand(200,200,10),[]);caxis([0 5]);
%Author:%Rohan, 1/3/2017, Matlab2016a.
%If there is no point in current plane initially it errors out - fix this.
hf=figure;
hf.Visible='on';
ha=axes('Parent',hf,'Tag','Axis');
ha.Units='normalized';
ha.Position=[0.1 0.1 0.8 0.8];
h_im=imagesc(IM(:,:,1),'Parent',ha);hold on

if ~isempty(r)
    hf.UserData.outplane.propname={'LineStyle','Marker','Color','LineWidth','MarkerSize'};
    hf.UserData.outplane.propval={'none','.',[0.5 0 0],0.5,8};
    hf.UserData.inplane.propname={'LineStyle','Marker','Color','LineWidth','MarkerSize'};
    hf.UserData.inplane.propval={'none','o',[1 0 0],1.5,8};
    h_points=gobjects(size(IM,3),1);
    for z=1:size(IM,3)
        incurrentplane=round(r(:,3))==z;
        h_points(z)=plot(r(incurrentplane,2),r(incurrentplane,1),'Tag',num2str(z),'Parent',ha);
        if z==1 %initial display has z=1 plane in view
            set(h_points(z),hf.UserData.inplane.propname,hf.UserData.inplane.propval);
        else
            set(h_points(z),hf.UserData.outplane.propname,hf.UserData.outplane.propval);
        end
    end
    hf.UserData.h_points=h_points;
end

hf.UserData.IM=IM;
hf.UserData.currplane=1;
hf.UserData.h_im=h_im;
hf.UserData.tt=tt;
ha.Title.String=['Current plane: ',num2str(hf.UserData.tt(hf.UserData.currplane)+1)];
hf.WindowScrollWheelFcn=@ScrollFcn;
%The figure has a property called WindowScrollWheelFcn that can be set to
%any function. A handle to such a function is passed, and the function is
%referred to as a callback function.
hf.Visible='on';
end

function ScrollFcn(src,ed)
%Callback function that is executed whenever the user scrolls. Callback
%functions that are executed when particular actions are performed (in this
%case a scroll operation) have 2 default inputs. Here the src is the figure
%handle, and ed is the event data. For WindowScrollWheelFcn, the default
%eventdata is a structure containing VerticalScrollCount along with other
%information.
hf=gcbf;
ha=findobj(hf.Children,'Type','Axes');
hi=hf.UserData.h_im;

%Update image
if ~isempty(hi)
    prevplane=hf.UserData.currplane;
    if ed.VerticalScrollCount>0
        hf.UserData.currplane=min(hf.UserData.currplane+ed.VerticalScrollCount,size(hf.UserData.IM,3));
    else
        hf.UserData.currplane=max(hf.UserData.currplane+ed.VerticalScrollCount,1);
    end
    hi.CData=hf.UserData.IM(:,:,hf.UserData.currplane);
end

%Update points
if isfield(hf.UserData,'h_points')
    set(hf.UserData.h_points(prevplane),hf.UserData.outplane.propname,hf.UserData.outplane.propval)
    set(hf.UserData.h_points(hf.UserData.currplane),hf.UserData.inplane.propname,hf.UserData.inplane.propval)
end

%Update axis
ha.Title.String=['Current plane: ',num2str(hf.UserData.tt(hf.UserData.currplane)+1)];
%ha.Title.String=['Current plane: ',num2str(hf.UserData.currplane)];
end