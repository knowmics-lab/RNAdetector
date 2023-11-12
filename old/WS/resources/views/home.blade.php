@extends('layouts.app')

@section('content')
    <div class="container">
        <div class="row justify-content-center">
            <div class="col-md-8">
                <div class="card">
                    <div class="card-header">Dashboard</div>
                    <div class="card-body">
                        @if (session('status'))
                            <div class="alert alert-success alert-dismissible fade show" role="alert">
                                {{ session('status') }}
                                <button type="button" class="close" data-dismiss="alert" aria-label="Close">
                                    <span aria-hidden="true">&times;</span>
                                </button>
                            </div>
                        @endif

                        <p class="card-text">Welcome <strong>{{ Auth::user()->name }}</strong>! Select an action:</p>

                        <div class="btn-toolbar justify-content-center" role="toolbar" aria-label="Actions toolbar">
                            <div class="btn-group mr-2" role="group" aria-label="Available actions">
                                @if (Auth::user()->admin)
                                    <a href="{{ route('users-list') }}" class="btn btn-secondary">Manage users</a>
                                    <a href="{{ route('run-update') }}" class="btn btn-secondary">Check db updates</a>
                                @endif
                                <a href="{{ route('change-password') }}" class="btn btn-secondary">Change Password</a>
                                <a href="{{ route('reset-token') }}" class="btn btn-secondary">Reset access token</a>
                                <a href="Javascript: void(0);" class="btn btn-secondary logout-button">Logout</a>
                            </div>
                        </div>

                    </div>
                </div>
            </div>
        </div>
        <div class="row justify-content-center mt-4">
            <div class="col-md-8">
                <div class="card">
                    <div class="card-header">
                        <div class="row">
                            <div class="col-3">
                                <i class="fas fa-5x fa-cogs"></i>
                            </div>
                            <div class="col-9 text-right">
                                <div class="huge">{{ $stats['jobs']['all'] }}</div>
                                <div>Jobs</div>
                            </div>
                        </div>
                    </div>
                    <ul class="list-group">
                        <li class="list-group-item">
                            <i class="fa fa-pause"></i>
                            Ready: <span class="badge badge-secondary badge-pill">{{ $stats['jobs']['ready'] }}</span>
                        </li>
                        <li class="list-group-item">
                            <i class="far fa-hourglass"></i>
                            Queued: <span class="badge badge-dark badge-pill">{{ $stats['jobs']['queued'] }}</span>
                        </li>
                        <li class="list-group-item">
                            <i class="fa fa-sync fa-spin"></i>
                            Processing: <span
                                class="badge badge-info badge-pill">{{ $stats['jobs']['processing'] }}</span>
                        </li>
                        <li class="list-group-item">
                            <i class="fa fa-check-circle"></i>
                            Completed: <span class="badge badge-danger badge-pill">{{ $stats['jobs']['completed'] }}</span>
                        </li>
                        <li class="list-group-item">
                            <i class="fa fa-exclamation-circle"></i>
                            Failed: <span
                                class="badge badge-success badge-pill">{{ $stats['jobs']['failed'] }}</span>
                        </li>
                    </ul>
                </div>
            </div>
        </div>
    </div>
@endsection
