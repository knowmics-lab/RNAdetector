<?php

namespace App\Http\Controllers;

use App\Models\User;
use Carbon\Carbon;
use Hash;
use Illuminate\Contracts\Support\Renderable;
use Illuminate\Http\RedirectResponse;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;

class UserController extends Controller
{

    /**
     * Create a new controller instance.
     *
     * @return void
     */
    public function __construct()
    {
        $this->middleware('auth');
    }

    /**
     * Shows the list of users
     *
     * @return \Illuminate\Contracts\View\Factory|\Illuminate\View\View
     */
    public function index(): Renderable
    {
        $users = User::paginate();

        return view(
            'admin.user.index',
            [
                'users' => $users,
            ]
        );
    }

    /**
     * Shows the user creation form
     *
     * @return \Illuminate\Contracts\Support\Renderable
     */
    public function create(): Renderable
    {
        return view('admin.user.create');
    }

    /**
     * Create a new user
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return \Illuminate\Http\RedirectResponse
     * @throws \Illuminate\Validation\ValidationException
     */
    public function doCreate(Request $request): RedirectResponse
    {
        $values = $this->validate(
            $request,
            [
                'name'     => ['required', 'max:255'],
                'email'    => ['required', 'max:255', 'email', 'unique:users'],
                'password' => ['required', 'string', 'min:6', 'alpha_num'],
                'admin'    => ['filled', 'boolean'],
            ]
        );
        $model = User::create(
            [
                'name'      => $values['name'],
                'email'     => $values['email'],
                'password'  => Hash::make($values['password']),
                'admin'     => $values['admin'] ?? false,
                'api_token' => null,
            ]
        );
        $model->email_verified_at = Carbon::now();
        $model->save();

        return redirect()->route('users-list')->with('status', 'User has been created.');
    }

    /**
     * Shows the user update form
     *
     * @param \App\Models\User $user
     *
     * @return \Illuminate\Contracts\Support\Renderable
     */
    public function show(User $user): Renderable
    {
        return view('admin.user.show', ['user' => $user]);
    }

    /**
     * Update the user
     *
     * @param \Illuminate\Http\Request $request
     * @param \App\Models\User         $user
     *
     * @return \Illuminate\Http\RedirectResponse
     * @throws \Illuminate\Validation\ValidationException
     */
    public function update(Request $request, User $user): RedirectResponse
    {
        $values = $this->validate(
            $request,
            [
                'name'     => ['required', 'max:255'],
                'email'    => ['required', 'max:255', 'email', Rule::unique('users')->ignore($user->id)],
                'password' => ['nullable', 'string', 'min:6', 'alpha_num'],
                'admin'    => ['filled', 'boolean'],
            ]
        );
        $user->name = $values['name'];
        $user->email = $values['email'];
        if (isset($values['password']) && !empty($values['password'])) {
            $user->password = Hash::make($values['password']);
        }
        $user->admin = (bool)($values['admin'] ?? false);
        $user->save();

        return redirect()->route('users-list')->with('status', 'User has been modified.');
    }

    /**
     * Delete the user
     *
     * @param \App\Models\User $user
     *
     * @return \Illuminate\Http\RedirectResponse
     * @throws \Exception
     */
    public function delete(User $user): RedirectResponse
    {
        $user->delete();

        return redirect()->route('users-list')->with('status', 'User has been deleted.');
    }

}
