function field = fieldFromComponents(xComponent,yComponent,zComponent)
%FIELDFROMCOMPONENTS Combine three scalar field interpolants into one vector field.
arguments
    xComponent (1,1) function_handle
    yComponent (1,1) function_handle
    zComponent (1,1) function_handle
end
field = @(x,y,z) [xComponent(x,y,z);yComponent(x,y,z);zComponent(x,y,z)];
end
