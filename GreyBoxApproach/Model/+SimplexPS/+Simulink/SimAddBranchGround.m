% This function addes branch ground.

% Author(s): Yitong Li

function SimAddBranchGround(Name_Model,Size_B_GND,Shift_B_GND,FullName_Branch,Name_Branch,ListLine)

fb = ListLine(:,1); % From bus
tb = ListLine(:,2); % To bus
N_Branch = length(fb);

for i = 1:N_Branch
    if fb(i) == tb(i)
        Name_B_GND{i} = ['B-GND' num2str(i)];
        FullName_BranchGND{i} = [Name_Model '/' Name_B_GND{i}];
        add_block('powerlib/Elements/Ground',FullName_BranchGND{i});
        PortPos_Branch{i} = get_param(FullName_Branch{i},'PortConnectivity');
        Pos_B_GND{i} = PortPos_Branch{i}(5).Position;
        Pos_B_GND{i} = Pos_B_GND{i} + Shift_B_GND;
        set_param(FullName_BranchGND{i},'position',[Pos_B_GND{i},Pos_B_GND{i}+Size_B_GND]);
        add_line(Name_Model,[Name_Branch{i} '/RConn2'],[Name_B_GND{i} '/LConn1'], ...
            'autorouting','smart');
    end
end

end