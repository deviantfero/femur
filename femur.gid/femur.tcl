proc InitGIDProject { dir } {
	Femur::SetDir $dir
	Femur::CreateWindow
}

# Defining namespace

namespace eval Femur {
	variable problemtype_dir;
}

proc Femur::SetDir { dir } {
	variable problemtype_dir
	set problemtype_dir $dir
}

proc Femur::GetDir { dir } {
	variable problemtype_dir
	return $problemtype_dir
}

proc Femur::CreateWindow { } {
	if { [GidUtils::AreWindowsDisabled] } {
		return
	}

	set w .gid.win_example
	InitWindow $w [= "FEMUR"] Femur "" "" 1

    if { ![winfo exists $w] } return ;# windows disabled || usemorewindows == 0
	ttk::frame $w.top
	ttk::label $w.top.title_text -text [= "FEMUR Window"]

	ttk::frame $w.bottom
	ttk::button $w.bottom.start -text [= "Start" ] -command [list destroy $w]

    grid $w.top.title_text -sticky ew
    grid $w.top -sticky new   

    grid $w.bottom.start -padx 6
    grid $w.bottom -sticky sew -padx 6 -pady 6
}
